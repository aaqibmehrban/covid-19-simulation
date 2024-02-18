import plotly.graph_objects as go
import numpy as np
from scipy import stats
from simulation import *
import copy

 # The flow matrix data
f = open('data/pij.pkl', 'rb')
flow_matrix_2020 = pickle.load(f)
f.close()

# The id of the city data
f = open('data/nodes.pkl', 'rb')
nodes_2020 = pickle.load(f)
f.close()

# The info of the city data
f = open('data/city_info.pkl', 'rb')
city_properties_2020 = pickle.load(f)
f.close()

first_cases = 0
firsts_day = []
firsts_dis = []
city_names_dis = []

# Get the virus info of different cities
shanghai_virus_info = preprocess_virus_info("上海市")
wuhan_virus_info = preprocess_virus_info("武汉市")



def comparison_date_distance_plotly():
    yy = np.array(firsts_day)
    xx = np.array(firsts_dis)

    bools_default = (xx != np.inf)
    bools = (xx != np.inf) & (xx < 5.5)
    slope, intercept, r_value, _, _ = stats.linregress(xx[bools], yy[bools])

    bools1 = (xx != np.inf) & (xx > 5.5)
    slope1, intercept1, r_value1, _, _ = stats.linregress(xx[bools1], yy[bools1])

    # Create figure
    fig = go.Figure()

    # Add traces
    fig.add_trace(go.Scatter(x=xx[bools], y=yy[bools], mode='markers', marker_color='royalblue', name='Effective Distance < 5.5'))
    fig.add_trace(go.Scatter(x=xx[bools1], y=yy[bools1], mode='markers', marker_color='orange', name='Effective Distance > 5.5'))
    fig.add_trace(go.Scatter(x=[firsts_dis[41]], y=[firsts_day[41]], mode='markers', marker_color='firebrick', name='Wuhan', marker_size=10))

    # Add lines
    fig.add_trace(go.Scatter(x=xx[bools_default], y=slope * xx[bools_default] + intercept, mode='lines', name='Linear Fit: All Cities', line=dict(color='firebrick')))
    fig.add_trace(go.Scatter(x=xx[bools1], y=slope1 * xx[bools1] + intercept1, mode='lines', name='Linear Fit: Distances > 5.5', line=dict(color='forestgreen')))

    # Annotations
    annotations = [
        dict(x=firsts_dis[0], y=firsts_day[0], xref='x', yref='y', text='Suzhou', showarrow=True, arrowhead=1, ax=0, ay=-40),
        dict(x=firsts_dis[11], y=firsts_day[11], xref='x', yref='y', text='Beijing', showarrow=True, arrowhead=1, ax=0, ay=-40),
        dict(x=firsts_dis[15], y=firsts_day[15], xref='x', yref='y', text='Chongqing', showarrow=True, arrowhead=1, ax=0, ay=-40),
        dict(x=firsts_dis[41], y=firsts_day[41], xref='x', yref='y', text='Wuhan', showarrow=True, arrowhead=1, ax=0, ay=-40)
    ]
    fig.update_layout(annotations=annotations)

    # Update layout
    fig.update_layout(title='Comparison of Reported Date After First Patient by Effective Distance',
                      xaxis_title='Effective Distance',
                      yaxis_title='Reported Date After First Patient',
                      legend_title='Legend',
                      legend=dict(orientation="h"),
                      margin=dict(l=0, r=0, t=30, b=0))

    return fig



def calculate_effective_distance_plotly(flow_matrix_data, nodes, virus_info, eff_distance=None):
    global firsts_day, firsts_dis
    if eff_distance is None:
        try:
            loaded_distance = np.load("data/effective_distance.npy")
            if loaded_distance.size == 0:
                effective_distance = eff_distance(flow_matrix_data)  # Assuming eff_distance() is defined elsewhere
            else:
                effective_distance = loaded_distance
        except FileNotFoundError:
            effective_distance = eff_distance(flow_matrix_data)  # Assuming eff_distance() is defined elsewhere
    else:
        effective_distance = eff_distance

    eff_dist = copy.deepcopy(effective_distance)
    for i in range(len(nodes)):
        eff_dist[i, i] = np.inf

    # Identify Shanghai's index
    shanghai_id = nodes['上海市']
    sorted_distance = np.sort(eff_dist[shanghai_id, :])
    sorted_distance_cities = np.argsort(eff_dist[shanghai_id, :])

    firsts_day = []
    firsts_dis = []
    city_names_dis = []

    for n, i in enumerate(sorted_distance_cities):
        city_name = find_in_nodes(i, nodes)  # Assuming find_in_nodes() is defined elsewhere
        item = virus_info.get(city_name, [])
        if len(item) == 0:
            item = virus_info.get(city_name[:-1], [])
        if len(item) > 0:
            firsts_day.append(item[0][0])
            firsts_dis.append(sorted_distance[n])
            city_names_dis.append(city_name)

    # Create Plotly figure
    fig = go.Figure()

    fig.add_trace(go.Scatter(x=sorted_distance, y=flow_matrix_data[shanghai_id, sorted_distance_cities],
                             mode='markers', marker=dict(color='royalblue', size=5),
                             name='Population Flow'))

    # Update layout for aesthetics similar to the original Matplotlib plot
    fig.update_layout(
        title='Effective Distance vs. Proportion of Population Flow',
        xaxis_title='Effective Distance',
        yaxis_title='Proportion of Population Flow',
        yaxis_type='log',  # Set y-axis to logarithmic scale
        xaxis=dict(range=[2.5, 15]),
        margin=dict(l=0, r=0, t=30, b=0),
        legend_title='Legend',
        legend=dict(orientation="h")
    )

    # Make top and right axis lines invisible
    fig.update_xaxes(showline=True, linewidth=2, linecolor='black', mirror=True)
    fig.update_yaxes(showline=True, linewidth=2, linecolor='black', mirror=True)

    return fig