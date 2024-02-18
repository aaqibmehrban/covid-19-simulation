from dash import Dash, html, dcc
import plotly.express as px
from simulation import *
from updated_graphs import*








fig1=calculate_effective_distance_plotly(flow_matrix_2020, nodes_2020, shanghai_virus_info, eff_distance=None)
fig2 = comparison_date_distance_plotly()
fig3=sir.apply_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, False)
fig4=sir.apply_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, True)
fig5=sir.verify_sir_model_plotly(flow_matrix_2020, nodes_2020, city_properties_2020, wuhan_virus_info, first_cases, False)


app = Dash(__name__)
server = app.server

app.layout = html.Div([
    html.Div([
        html.H1("Covid 19 Simulation", style={'textAlign': 'center'})
    ], style={'background': 'linear-gradient(to right, lightgray, white)', 'padding': '10px', 'margin-bottom': '20px'}),
    
    # Body
    html.Div([
        # Graphs container
        html.Div([
            # Column 1
            html.Div([
                html.Div([dcc.Graph(id='graph-1', figure=fig1)], style={'background-color': '#f8f9fa', 'margin': '10px', 'padding': '20px', 'border-radius': '5px'}),
                html.Div([dcc.Graph(id='graph-3', figure=fig3)], style={'background-color': '#f8f9fa', 'margin': '10px', 'padding': '20px', 'border-radius': '5px'}),
                html.Div([dcc.Graph(id='graph-5', figure=fig5)], style={'background-color': '#f8f9fa', 'margin': '10px', 'padding': '20px', 'border-radius': '5px'})
            ], style={'width': '50%', 'display': 'inline-block', 'verticalAlign': 'top'}),

            # Column 2
            html.Div([
                html.Div([dcc.Graph(id='graph-2', figure=fig2)], style={'background-color': '#f8f9fa', 'margin': '10px', 'padding': '20px', 'border-radius': '5px'}),
                html.Div([dcc.Graph(id='graph-4', figure=fig4)], style={'background-color': '#f8f9fa', 'margin': '10px', 'padding': '20px', 'border-radius': '5px'})
            ], style={'width': '50%', 'display': 'inline-block', 'verticalAlign': 'top'})
        ], style={'display': 'flex', 'flexDirection': 'row'})
    ], style={'padding': '20px'})
], style={'background-color': 'white', 'fontFamily': 'Arial, sans-serif'})




if __name__ == '__main__':
    port = int(os.environ.get('PORT', 5000))
    app.run_server(host='0.0.0.0', port=port)
    

