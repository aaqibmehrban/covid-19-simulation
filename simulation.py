import math
import copy
import pickle
import sir_model as sir
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from tqdm import trange


# Get the virus info/infected numbers from the news
# As the data is mainly in Chinese, and we do not have so much time to translate it
# So it is sorry to make it not readable
# 上海市 - Shanghai
# 武汉市 - Wuhan
def preprocess_virus_info(city):
    # Load the data
    news_data = pd.ExcelFile('data/city_cases.xlsx').parse("Sheet1")
    # Retrieve unique cities from the data, excluding the first entry
    all_reported_cities = list(set(news_data['城市']))[1:]
    # Extend the city list with additional key cities
    all_reported_cities = all_reported_cities + ['北京市', '深圳市', '重庆市', '广州市']
    # Filter data for the specified city
    city_virus = news_data.loc[
        news_data['城市'] == city, ['新增确诊病例', '确诊/出院', '公开报道时间', '新增治愈出院数', '新增死亡数']]

    # Sort the data by diagnosis/outcome date
    city_virus.sort_values(by='确诊/出院')
    # Reset index after sorting
    city_virus.reset_index(drop=True, inplace=True)
    # Get the date of the first reported case
    first_date = city_virus.iloc[-1]['确诊/出院']
    # Initialize global variable for the first cases count
    global first_cases
    # Determine the initial number of cases, using a minimum of 4
    first_cases = 4 if 4 > city_virus.iloc[-1]['新增确诊病例'] else city_virus.iloc[-1]['新增确诊病例']

    # municipality is used to distinguish the municipality
    # which means the city name is the province name
    all_cities_cases = {}
    for city in all_reported_cities:
        subset = news_data.loc[
            news_data['城市'] == city, ['新增确诊病例', '确诊/出院', '公开报道时间', '新增治愈出院数', '新增死亡数',
                                        '省份']]
        municipality = False
        if len(subset) == 0:
            subset = news_data.loc[
                news_data['省份'] == city[:-1], ['新增确诊病例', '确诊/出院', '公开报道时间', '新增治愈出院数',
                                                 '新增死亡数', '省份']]
            municipality = True

        # Extract arrays of cases, cures, deaths and dates
        new_cases = np.array(subset["新增确诊病例"])
        cued_cases = np.array(subset['新增治愈出院数'])
        die_cases = np.array(subset['新增死亡数'])
        found_virus_dates = list(subset['确诊/出院'])
        reported_dates = list(subset['公开报道时间'])
        # Calculate days from the first reported case
        days = []
        for i, date in enumerate(found_virus_dates):
            # Use reported date if diagnosis date is missing
            if pd.isna(date):
                date = reported_dates[i]
            if not pd.isna(date):
                days.append(int((date - first_date) / np.timedelta64(1, 'D')))

        sorted_days = np.sort(days)
        # Statistic of patients in each city during the period(2020.1.10 - 2020.1.29)
        days_index = np.argsort(days)
        infected = np.cumsum(new_cases[days_index])
        cued = np.cumsum(cued_cases[days_index])
        death = np.cumsum(die_cases[days_index])
        if len(sorted_days) > 0:
            if municipality:
                all_cities_cases[list(subset['省份'])[0] + '市'] = (sorted_days, infected, cued, death)
            else:
                all_cities_cases[city] = (sorted_days, infected, cued, death)
    return all_cities_cases


def my_matrix_production(matrix, matrix_copy):
    # Perform a special kind of matrix multiplication that uses 'maximum' instead of 'sum' across the third dimension.
    # Size of the matrix (assuming square matrix)
    sz = matrix.shape[0]
    # Initialize the output matrix with zeros
    output = np.zeros([sz, sz])
    # Perform the custom matrix product
    for i in range(sz):
        output[i, :] = np.max(matrix[i, :, None] * matrix_copy, axis=0)
    return output


def eff_distance(prob):
    # Calculate the effective distance between nodes in a network based on a probability matrix.
    # Size of the probability matrix (assuming square matrix)
    sz = prob.shape[0]
    # Small constant to prevent log of zero
    epsilon = 1e-9
    # Start with the original probability matrix
    prod = prob
    # Initialize the distance matrix based on the probability matrix
    distance = np.ones([sz, sz]) - np.log(prod + epsilon)
    # Iteratively calculate the effective distance
    for i in trange(1, sz - 1):
        # Update the product matrix using the custom matrix product function
        prod = my_matrix_production(prod, prob)
        # Calculate the distance matrix for this iteration
        dist = i + 1 - np.log(prod + epsilon)
        # Update the overall distance matrix with the minimum distances
        distance = np.minimum(distance, dist)
    # Save the effective distance matrix to a file
    np.save("data/effective_distance.npy", distance)
    return distance


def calculate_effective_distance(flow_matrix_data, nodes, virus_info):
    try:
        loaded_distance = np.load("data/effective_distance.npy")
        if loaded_distance.size == 0:
            effective_distance = eff_distance(flow_matrix_data)
        else:
            effective_distance = loaded_distance
    except FileNotFoundError:
        effective_distance = eff_distance(flow_matrix_data)
    eff_dist = copy.deepcopy(effective_distance)
    for i in range(len(nodes)):
        eff_dist[i, i] = np.inf

    # Sort by effective distance
    shanghai_id = nodes['上海市']
    sorted_distance = np.sort(eff_dist[shanghai_id, :])
    sorted_distance_cities = np.argsort(eff_dist[shanghai_id, :])

    # Sort by flow data
    # sorted_flux = np.sort(flow_matrix_data[shanghai_id, :])[::-1]
    # sorted_flux_cities = np.argsort(flow_matrix_data[shanghai_id, :])[::-1]

    for n, i in enumerate(sorted_distance_cities):
        city_name = find_in_nodes(i, nodes)
        item = virus_info.get(city_name, [])
        # item : sorted_days, infected, cued, death
        if len(item) == 0:
            item = virus_info.get(city_name[:-1], [])
        if len(item) > 0:
            firsts_day.append(item[0][0])
            firsts_dis.append(sorted_distance[n])
            city_names_dis.append(city_name)
        # print(city_name, sorted_distance[n])

    plt.figure(figsize=(10, 6))
    # This turns off minor ticks on the y-axis
    plt.tick_params(axis='y', which='minor', left=False)
    plt.xlim(2.5, 15)
    plt.rcParams['axes.unicode_minus'] = False
    plt.plot(eff_dist[shanghai_id, :], flow_matrix_data[shanghai_id, :], 'o', color='royalblue', markersize=5)
    plt.xlabel('Effective Distance')
    plt.ylabel('Proportion of Population Flow')
    plt.yscale('log')
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('figure/effective_distance_flow.png', dpi=600)
    plt.show()


def find_in_nodes(index, nodes_set):
    key = -1
    for k, value in nodes_set.items():
        if value == index:
            key = k
            break
    return key


def comparison_date_distance():
    yy = np.array(firsts_day)
    xx = np.array(firsts_dis)
    plt.figure(figsize=(15, 10))

    bools_default = (xx != math.inf)
    bools = (xx != math.inf) & (xx < 5.5)
    slope, intercept, r_value, _, _ = stats.linregress(xx[bools], yy[bools])

    bools1 = (xx != math.inf) & (xx > 5.5)
    slope1, intercept1, r_value1, _, _ = stats.linregress(xx[bools1], yy[bools1])

    plt.scatter(xx[bools], yy[bools], s=60, color='royalblue')
    plt.scatter(xx[bools1], yy[bools1], s=12, color='orange')
    plt.scatter(firsts_dis[41], firsts_day[41], s=80, color='firebrick')

    plt.plot(xx[bools_default], slope * xx[bools_default] + intercept, color='firebrick', label='All the cities')
    plt.plot(xx[bools1], slope1 * xx[bools1] + intercept1, color='forestgreen',
             label='Only cities with effective distances > 5.5')

    # 0 Suzhou 8 Yancheng 11 Beijing 15 Chongqing
    plt.annotate("Suzhou", (firsts_dis[0], firsts_day[0]), textcoords="offset points", xytext=(0, 10),
                 ha='center', fontsize=14)
    plt.annotate("Beijing", (firsts_dis[11], firsts_day[11]), textcoords="offset points", xytext=(0, 10),
                 ha='center', fontsize=14)
    plt.annotate("Chongqing", (firsts_dis[15], firsts_day[15]), textcoords="offset points", xytext=(0, 10),
                 ha='center', fontsize=14)
    plt.annotate("Wuhan", (firsts_dis[41], firsts_day[41]), textcoords="offset points", xytext=(0, 10),
                 ha='center', fontsize=14)

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.rcParams['axes.unicode_minus'] = False
    plt.xlabel('Effective Distance', fontsize=14)
    plt.ylabel('Reported Date After First Patient', fontsize=14)
    plt.legend(loc='upper left', shadow=True, numpoints=1, fontsize=10)
    plt.savefig('figure/effective_distance_date.png', dpi=600)
    plt.show()


if __name__ == '__main__':
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

    # Calculate the effective distance
    calculate_effective_distance(flow_matrix_2020, nodes_2020, shanghai_virus_info)

    # Combine the date with the effective distance
    comparison_date_distance()

    # Apply the sir model for Shanghai
    sir.apply_sir_model(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, False)
    sir.apply_sir_model(flow_matrix_2020, nodes_2020, city_properties_2020, first_cases, True)
    sir.verify_sir_model(flow_matrix_2020, nodes_2020, city_properties_2020, wuhan_virus_info, first_cases, False)
