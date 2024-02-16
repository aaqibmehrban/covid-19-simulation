import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def protect_decay(t, t0, eta, rate_time):
    # Small constant to avoid division by zero
    epsilon = 0.001
    # Calculate decay rate
    r = 2 * np.log((1 - epsilon) / epsilon) / rate_time
    # Midpoint of decay
    x0 = t0 + rate_time / 2
    # Decay function for protective measures
    decay = eta / (1 + np.exp(r * (t - x0))) + 1 - eta
    return decay


def sir_model(population, t, alpha, beta, gamma, eta, rate_time, protect_day, flow_matrix_data, intervention):
    # Size of the susceptible and infected populations
    sz = population.shape[0] // 2
    # Infected population
    js = population[:sz]
    # Susceptible population
    ss = population[sz:]
    # Terms for movement of infected and susceptible populations
    jterm = js.dot(flow_matrix_data) - js
    sterm = ss.dot(flow_matrix_data) - ss
    # Cross term for infection rate, modified by protective measures if intervention is True
    if not intervention:
        cross_term = alpha * js * ss
    else:
        cross_term = alpha * protect_decay(t, protect_day, eta, rate_time) * js * ss
    # Change in infected population
    delta_i = cross_term - beta * js + gamma * jterm
    # Change in susceptible population
    delta_s = - cross_term + gamma * sterm
    # Combined change for the SIR model
    result = np.r_[delta_i, delta_s]
    return result


def verify_sir_model(flow_matrix_data, nodes, city_properties, virus_info, first_cases, intervention):
    timespan = np.linspace(0, 40, 80)
    js0 = np.zeros(len(nodes))
    ss0 = np.ones(len(nodes))
    city_initial = '武汉市'
    js0[nodes[city_initial]] = float(first_cases) / float(city_properties[city_initial]['pop'])
    r0 = 3.8
    beta = 1 / 10  # The time for the recovery
    alpha = r0 * beta  # Infection  rate
    gamma = 0.1255  # Flow rate
    eta = (1 - 4.0 * beta / r0) if intervention else (1 - 1.0 / r0)
    rate_time = 30  # Adjust rate

    # Apply SIR model
    simulation_result = odeint(sir_model, np.r_[js0, ss0], timespan,
                               args=(alpha, beta, gamma, eta, rate_time, 13,
                                     np.transpose(flow_matrix_data), intervention))

    cities = ['武汉市', '重庆市', '深圳市', '上海市', '北京市', '广州市']
    cities_name_eng = {'武汉市': 'Wuhan', '重庆市': 'Chongqing', '深圳市': 'Shenzhen', '北京市': 'Beijing',
                       '上海市': 'Shanghai', '广州市': 'Guangzhou'}

    plt.figure(figsize=(15, 8))
    plot_time_span = 700
    line_width_value = 2
    colors = plt.cm.jet(np.linspace(0, 1, len(cities)))

    # Plot city-specific SIR Model results and actual infection data
    for n, city_name in enumerate(cities):
        idx = nodes.get(city_name, -1)
        if idx >= 0:
            plt.plot(timespan[:plot_time_span],
                     simulation_result[:plot_time_span, idx] * city_properties[city_name]['pop'], color=colors[n],
                     label=cities_name_eng[city_name], linewidth=line_width_value)
            item = virus_info.get(city_name)
            if item:
                # infected_number = infected - cued - death
                infected_number = item[1] - item[2] - item[3]
                plt.plot(item[0], infected_number, 'o', color=colors[n])

    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend(loc='upper left', shadow=True, numpoints=1, fontsize=10)
    plt.xlabel('Days(From 2020.1.10)')
    plt.ylabel('Infected Population(Log)')
    plt.yscale('log')
    plt.savefig('figure/sir_model_verification_intervention_{}.png'.format(intervention), dpi=600)
    plt.show()


def apply_sir_model(flow_matrix_data, nodes, city_properties, first_cases, intervention):
    timespan = np.linspace(0, 200, 1000)
    js0 = np.zeros(len(nodes))
    ss0 = np.ones(len(nodes))
    city_initial = '上海市'
    js0[nodes[city_initial]] = float(first_cases) / float(city_properties[city_initial]['pop'])
    r0 = 3.8
    beta = 1 / 10  # The time for the recovery
    alpha = r0 * beta  # Infection  rate
    gamma = 0.1255  # Flow rate
    eta = (1 - 4.0 * beta / r0) if intervention else (1 - 1.0 / r0)
    rate_time = 30  # Adjust rate

    # Apply SIR model
    simulation_result = odeint(sir_model, np.r_[js0, ss0], timespan,
                               args=(alpha, beta, gamma, eta, rate_time, 13,
                                     np.transpose(flow_matrix_data), intervention))
    cities = ['上海市', '苏州市', '重庆市', '北京市', '盐城市', '南通市', '南京市', '广州市']
    cities_name_eng = {'上海市': 'Shanghai', '苏州市': 'Suzhou', '重庆市': 'Chongqing', '北京市': 'Beijing',
                       '盐城市': 'Yancheng', '南通市': 'Nantong', '南京市': 'Nanjing', '广州市': 'Guangzhou'}
    plt.figure(figsize=(15, 8))
    plot_time_span = 700
    alpha_value = 0.1
    line_width_value = 2
    colors = plt.cm.jet(np.linspace(0, 1, len(cities)))

    # Plot SIR Model results for all cities
    for i, city_name in enumerate(nodes):
        population = city_properties[city_name]['pop']
        plt.plot(timespan[:plot_time_span], simulation_result[:plot_time_span, i] * population, alpha=alpha_value)

    # Plot city-specific SIR Model results and actual infection data
    for n, city_name in enumerate(cities):
        idx = nodes.get(city_name, -1)
        if idx >= 0:
            plt.plot(timespan[:plot_time_span],
                     simulation_result[:plot_time_span, idx] * city_properties[city_name]['pop'], color=colors[n],
                     label=cities_name_eng[city_name], linewidth=line_width_value)
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.legend(loc='upper left', shadow=True, numpoints=1, fontsize=10)
    plt.xlabel('Days(From 2020.1.20)')
    plt.ylabel('Infected Population')
    plt.savefig('figure/sir_model_intervention_{}.png'.format(intervention), dpi=600)
    plt.show()
