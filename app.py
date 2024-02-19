import os

from dash import Dash, html, dcc
from updated_graphs import*



fig1=calculate_effective_distance_plotly(flow_matrix_2020, nodes_2020, shanghai_virus_info, eff_distance=None)
fig2 = comparison_date_distance_plotly()
fig3=sir.apply_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, False)
fig4=sir.apply_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, True)
fig5=sir.verify_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, wuhan_virus_info, first_cases, False)


app = Dash(__name__)
server=app.server


app.layout = html.Div([
    html.Div([
        html.H1("Covid 19 Simulation", style={'textAlign': 'center'})
    ], style={'background': 'linear-gradient(to right, lightgray, white)', 'padding': '1vh', 'margin-bottom': '2vh'}),

    html.Div([
        dcc.Graph(id='graph-1', figure=fig1, style={'width': '90vw', 'height': '45vh'}),
        dcc.Graph(id='graph-2', figure=fig2, style={'width': '90vw', 'height': '45vh'}),
        dcc.Graph(id='graph-3', figure=fig3, style={'width': '90vw', 'height': '45vh'}),
        dcc.Graph(id='graph-4', figure=fig4, style={'width': '90vw', 'height': '45vh'}),
        dcc.Graph(id='graph-5', figure=fig5, style={'width': '90vw', 'height': '45vh'}),
    ], style={'textAlign': 'center'})
], style={'background-color': 'white', 'fontFamily': 'Arial, sans-serif', 'maxWidth': '100vw'})


if __name__ == '__main__':
    port = int(os.environ.get("PORT", 5000))
    app.run_server(host='0.0.0.0', port=port)

