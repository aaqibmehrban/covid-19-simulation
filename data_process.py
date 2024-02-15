import pickle
import pandas as pd


# Read the city id info and store it into the nodes
def id_to_city(csv_file):
    city_info = {}
    for i in csv_file.iterrows():
        city_info[i[1]['CITY_ID']] = {'city': i[1]['CITY'], 'prov': i[1]['PROV']}
    return city_info


# Transform the flow data into the form of matrix
def read_flow_data():
    flow_data = pd.read_csv('Pij_BAIDU.csv', encoding='gbk')
    # slices of the method will automatically avoid the first row
    city_flow_matrix = flow_data.iloc[:, 1:].values / 100
    return city_flow_matrix


# The whole workflow of the data preprocess
# The process has abandoned as we found the data which has already been processed
def preprocess_data_workflow():
    city_info_2020 = id_to_city(pd.read_csv('county_city_province.csv'))
    city_flow_matrix_2020 = read_flow_data()
    pickle.dump(city_info_2020, open('city_info_2020.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)
    pickle.dump(city_flow_matrix_2020, open('city_flow_matrix_2020.pkl', 'wb'), pickle.HIGHEST_PROTOCOL)