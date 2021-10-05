# this template have two connections including one Power Grid connection
# len(g_connection[i]) == 2
# gas system: one well(Grid) one well(Connection)
import numpy as np

a = 2

T = 2             # this is power time slot
time_point = 500
delta_time = 1
node_voltage_min = 0.6
node_voltage_max = 1.6
gas_flow_min = 0
gas_flow_max = 10
gas_buy_price = 1

tc = 0.333            # 1 / tc = tao , tao * time_point = 300s, total = 300 * 10 = 3000 s
pipeline_length = 9000
w_av = 2.55
c1 = 300
diameter = 0.5
Lambda = 0.07 * 1
Time = T * time_point

power_load_num = 8

load_power_min_raw = np.array([1.69, 1.729, 1.664, 1.703, 1.742, 1.833, 1.807, 1.755, 1.768, 1.729])
load_power_max_raw = load_power_min_raw * 1.4
load_share = [1 / power_load_num] * power_load_num
load_power_min = np.array(load_share).reshape((-1, 1)).dot(load_power_min_raw.reshape((1, -1)))
load_power_min = load_power_min[:, 0:T].tolist()
load_power_max = np.array(load_share).reshape((-1, 1)).dot(load_power_max_raw.reshape((1, -1)))
load_power_max = load_power_max[:, 0:T].tolist()


##
#  (GRID)                                                                  (GRID)
# d-12-g o                                                                   0-w o
#        |↓                                                                      ↘╲
#    d-2 |↓ <-  1     <-    0-gg ->   d-3  ->    d-4             4-l      2  <-   ↘╲      ->     3      5-gas-gen
#        o--------o---------o----------o-----------o              o-------o-----------o----------o--------o
#       ▲                  |                                              ▼            1          ▼
#       d-7 <-    6   <-    |5-g ->  d-8   ->     d-9
#        o--------o---------o----------o-----------o
#               ↓ |       ↓ |
#               ↓ |       ↓ |
#           d-11  o         o d-10
# #

# change 1
# 输出功率的上下限
gen_power_min = [0,         0,         0,   ]
gen_power_max = [0.3,       0.1,       3.1, ]
gen_react_min = [0,         0,         0,   ]
gen_react_max = [0.3,       0.1,       3.1, ]
# 发电机成本
gen_cost_a = [   12.8875,   12.8875,   26.2438,    ]
gen_cost_b = [   0.010875,  0.010875,  0.0699663,  ]
gen_cost_c = [   1,         1,         1,          ]
gen_index =  [   0,         0,         5,          ]  # this is start

load_index = [   2, 3, 4, 7, 8, 9, 10, 11,         ]  # this is start
# load_index = [   2, 3]  # this is start
assert len(load_index) == power_load_num
upper_grid_connection_index_power = [1, ]  # this is end
upper_well_connection_index_gas =   [1, ]  # this is end

line_start_point = [0, 1, 0, 3, 0, 5, 8, 5, 6, 5,  6,  12, ]  # this is start ->
line_end_point =   [1, 2, 3, 4, 5, 8, 9, 6, 7, 10, 11, 1,  ]  # end
# line_start_point = [0, 1, 0, 3, 0, ]  # this is start ->
# line_end_point =   [1, 2, 3, 4, 5, ]  # end
line_current_flow_capacity = 10

#######################################################################
gas_node_num =              6
node_pressure_min =         4e4 / 4e3 * 5    # 50
node_pressure_max =         8e5 / 4e3        # 200
# node_pressure_max = 5e1 / 4e3
PRESSURE_CONSTANT =         4e5 / 4e3       # 100
# PRESSURE_CONSTANT = 4e0 / 4e3

gas_well_num =              1
well_min =                 [0, ]
well_max =                 [20, ]
well_index =               [0, ]  # 6


gas_load_index =           [4 ]  # 5
gas_load_num =             len(gas_load_index)



load_gas_min_raw = [(np.array([0.1, 0.21, 0.22, 0.11, 0.12, 0.11, 0.12, 0.11, 0.10, 0.10])*10).tolist()]
load_gas_max_raw = [(np.array([0.2, 0.41, 0.42, 0.21, 0.23, 0.24, 0.21, 0.26, 0.23, 0.24])*10).tolist()]
load_gas_min = (np.array(load_gas_min_raw))[:, 0:T].tolist() * gas_load_num
load_gas_max = (np.array(load_gas_max_raw))[:, 0:T].tolist() * gas_load_num




gas_line_num =             gas_node_num - 1
gas_line_start_point =     [0, 1, 1, 2, 3 ]
gas_line_end_point =       [1, 2, 3, 4, 5 ]

gen_gas_capacity =         [20, ]
gen_gas_power_min =        [0, ]

gen_gas_power_index =      [2, ]
gen_gas_gas_index =        [5, ]

gen_gas_efficiency =       [1, ]








line_resistance = (np.array(
    [0.1927, 0.0559, 0.0385, 0.05593, 0.0895, 0.1927, 0.03854, 0.03854, 0.1927, 0.0895, 0.0559, 0.0559  # local
     ])*0).tolist()
line_reactance = (np.array(
    [0.06562, 0.03563, 0.0376, 0.03563, 0.05701, 0.06562, 0.0376, 0.0376, 0.06562, 0.05701, 0.05701, 0.0559  # local
     ])*0).tolist()


system_info = {
    'T': T,
    'time_point': time_point,
    'delta_time': delta_time,
}

thermal_gen_info = {
    'gen_num': len(gen_power_min),
    'gen_index': gen_index,
    'gen_power_min': gen_power_min,
    'gen_power_max': gen_power_max,
    'gen_react_min': gen_react_min,
    'gen_react_max': gen_react_max,
    'gen_cost_a': gen_cost_a,
    'gen_cost_b': gen_cost_b,
    'gen_cost_c': gen_cost_c,

    'connection_index_power': upper_grid_connection_index_power,
    'connection_index_gas': upper_well_connection_index_gas,
}

bus_info = {
    'bus_num': len(line_start_point) + 1,
    'bus_voltage_min': [node_voltage_min] * (len(line_start_point) + 1),
    'bus_voltage_max': [node_voltage_max] * (len(line_start_point) + 1),
}

load_info = {
    'load_num': len(load_index),
    'load_index': load_index,
    'load_power_min': load_power_min,  #
    'load_power_max': load_power_max,  #
    'load_react_min': load_power_min,  #
    'load_react_max': load_power_max,  #
}

line_info = {
    'line_num': len(line_start_point),
    'line_current_capacity': [line_current_flow_capacity] * len(line_start_point),
    'line_power_flow_capacity': line_current_flow_capacity,
    'line_react_flow_capacity': line_current_flow_capacity,
    'line_start_point': line_start_point,
    'line_end_point': line_end_point,
    'line_resistance': np.array(line_resistance),
    'line_reactance': np.array(line_reactance)
}

gas_node_info = {
    'gas_node_num': gas_node_num,
    'gas_node_pressure_min': node_pressure_min,
    'gas_node_pressure_max': node_pressure_max,
}

gas_line_info = {
    'gas_line_num': gas_line_num,
    'gas_pipeline_start_node': gas_line_start_point,
    'gas_pipeline_end_node': gas_line_end_point,
    'gas_flow_max': gas_flow_max,
    'gas_flow_min': gas_flow_min,
}

gas_load_info = {
    'gas_load_num': gas_load_num,
    'gas_load_connection_index': gas_load_index,
    'gas_load_min': load_gas_min,
    'gas_load_max': load_gas_max,
}

gas_generator_info = {
    'gen_gas_num': len(gen_gas_capacity),
    'gen_gas_index_gas': gen_gas_gas_index,  # the gas generator index in the gas system
    'gen_gas_index_power': gen_gas_power_index,  # the gas generator index in the power system
    'gen_gas_min': gen_gas_power_min,  # this is power
    'gen_gas_capacity': gen_gas_capacity,
    'gen_gas_efficiency': gen_gas_efficiency,
}

gas_well_info = {
    'gas_well_num': gas_well_num,
    'gas_well_connection_index': well_index,
    'well_output_min': well_min,
    'well_output_max': well_max,
    'well_cost': gas_buy_price
}


def getConfig():
    power_config = {
        'thermal_gen_info': thermal_gen_info,
        'bus_info': bus_info,
        'power_load_info': load_info,
        'power_line_info': line_info
    }
    gas_config = {
        'gas_line_info' :gas_line_info,
        'gas_node_info' :gas_node_info,
        'gas_load_info' :gas_load_info,
        'gas_generator_info': gas_generator_info,
        'gas_well_info': gas_well_info,
        'PRESSURE_CONSTANT': PRESSURE_CONSTANT
    }
    pipeline_config = {
        'T_long' : T,
        'tc'     : tc,
        'L'      : pipeline_length,
        'w_av'   : w_av,
        'c1'     : c1,
        'dij'    : diameter,
        'Lambda' : Lambda,
        'Time_all'   : Time,
    }

    return power_config, gas_config, pipeline_config
















