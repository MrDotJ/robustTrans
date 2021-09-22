import numpy as np
import gurobipy as gurobi
from config.powergas import getConfig


def tonp(vars_gur) -> np.array:
    keys = vars_gur.keys()
    key_last = keys[-1]
    if type(key_last) != int:
        dim = len(key_last)
        shape = list(map(lambda x: x + 1, list(key_last)))
    else:
        dim = 1
        shape = key_last + 1
    if dim == 1:
        re = np.empty(tuple([shape, ]), dtype=object)
        for i in range(shape):
            re[i] = vars_gur[i]
        return re
    if dim == 2:
        re = np.empty(tuple(shape), dtype=object)
        for i in range(shape[0]):
            for j in range(shape[1]):
                re[i, j] = vars_gur[i, j]
        return re
    if dim == 3:
        re = np.empty(shape, dtype=object)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    re[i, j, k] = vars_gur[i, j, k]
        return re
    if dim == 4:
        re = np.empty(shape, dtype=object)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        re[i, j, k, m] = vars_gur[i, j, k, m]
        return re
    if dim == 5:
        re = np.empty(shape, dtype=object)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        for n in range(shape[4]):
                            re[i, j, k, m, n] = vars_gur[i, j, k, m, n]
        return re
def to_value(vars_gur):
    shape = vars_gur.shape
    dim = vars_gur.ndim
    if dim == 1:
        re = np.empty(tuple([shape[0], ]), dtype=np.float)
        for i in range(shape[0]):
            re[i] = vars_gur[i].getAttr('X')
        return re
    if dim == 2:
        re = np.empty(tuple(shape), dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                re[i, j] = vars_gur[i, j].getAttr('X')
        return re
    if dim == 3:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    re[i, j, k] = vars_gur[i, j, k].getAttr('X')
        return re
    if dim == 4:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        re[i, j, k, m] = vars_gur[i, j, k, m].getAttr('X')
        return re
    if dim == 5:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        for n in range(shape[4]):
                            re[i, j, k, m, n] = vars_gur[i, j, k, m, n].getAttr('X')
        return re
def to_value_tuple(vars_gur):
    keys = vars_gur.shape
    key_last = keys[-1]
    if type(key_last) != int:
        dim = len(key_last)
        shape = list(map(lambda x: x + 1, list(key_last)))
    else:
        dim = 1
        shape = key_last + 1
    if dim == 1:
        re = np.empty(tuple([shape, ]), dtype=np.float)
        for i in range(shape):
            re[i] = vars_gur[i].getAttr('X')
        return re
    if dim == 2:
        re = np.empty(tuple(shape), dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                re[i, j] = vars_gur[i, j].getAttr('X')
        return re
    if dim == 3:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    re[i, j, k] = vars_gur[i, j, k].getAttr('X')
        return re
    if dim == 4:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        re[i, j, k, m] = vars_gur[i, j, k, m].getAttr('X')
        return re
    if dim == 5:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    for m in range(shape[3]):
                        for n in range(shape[4]):
                            re[i, j, k, m, n] = vars_gur[i, j, k, m, n].getAttr('X')
        return re


class PowerGas:
    def __init__(self, powerConfig, gasConfig, pipelineConfig):
        pc = powerConfig
        gc = gasConfig
        thermal_gen_info = pc['thermal_gen_info']
        power_bus_info = pc['bus_info']
        power_load_info = pc['power_load_info']
        power_line_info = pc['power_line_info']
        gas_line_info = gc['gas_line_info']
        gas_node_info = gc['gas_node_info']
        gas_load_info = gc['gas_load_info']
        gas_generator_info = gc['gas_generator_info']
        gas_well_info = gc['gas_well_info']
        # ------------- generator of non-gas [power system]----------------------
        self.gen_num = thermal_gen_info['gen_num']  # add virtual node at last as the connected node
        self.gen_conn_power_index = np.array(thermal_gen_info['gen_index'])
        self.gen_power_min = thermal_gen_info['gen_power_min']
        self.gen_power_max = thermal_gen_info['gen_power_max']
        self.gen_react_min = thermal_gen_info['gen_react_min']
        self.gen_react_max = thermal_gen_info['gen_react_max']
        self.gen_cost_a = thermal_gen_info['gen_cost_a']
        self.gen_cost_b = thermal_gen_info['gen_cost_b']
        self.gen_cost_c = thermal_gen_info['gen_cost_c']
        # ---------------- power bus node ---------------------------
        self.bus_num = power_bus_info['bus_num']
        self.bus_voltage_min = power_bus_info['bus_voltage_min']
        self.bus_voltage_max = power_bus_info['bus_voltage_max']
        # ------------------- power line info -------------------------
        self.line_num = power_line_info['line_num']
        self.line_current_capacity = power_line_info['line_current_capacity']
        self.line_start_point = np.array(power_line_info['line_start_point'])
        self.line_end_point = np.array(power_line_info['line_end_point'])
        self.line_resistance = power_line_info['line_resistance']
        self.line_reactance = power_line_info['line_reactance']
        self.line_power_flow_max = power_line_info['line_power_flow_capacity']
        self.line_react_flow_max = power_line_info['line_react_flow_capacity']
        # ----------------- power load -------------------------------
        self.load_num = power_load_info['load_num']
        self.load_index = np.array(power_load_info['load_index'])
        self.load_power_min = np.array(power_load_info['load_power_min'])
        self.load_power_max = np.array(power_load_info['load_power_max'])
        self.load_react_min = np.array(power_load_info['load_react_min'])
        self.load_react_max = np.array(power_load_info['load_react_max'])
        # ----------------- gas line info [gas system]------------------------------
        self.gas_line_num = gas_line_info['gas_line_num']
        self.gas_flow_min = gas_line_info['gas_flow_min']
        self.gas_flow_max = gas_line_info['gas_flow_max']
        self.gas_pipeline_start_node = np.array(gas_line_info['gas_pipeline_start_node'])
        self.gas_pipeline_end_node = np.array(gas_line_info['gas_pipeline_end_node'])
        # ----------------- gas node info ----------------------------
        self.PRESSURE_CONSTANT = gc['PRESSURE_CONSTANT']
        self.gas_node_num = gas_node_info['gas_node_num']
        self.gas_node_pressure_min = gas_node_info['gas_node_pressure_min']
        self.gas_node_pressure_max = gas_node_info['gas_node_pressure_max']
        # ----------------- gas load info -----------------------------
        self.gas_load_num = gas_load_info['gas_load_num']
        self.gas_load_min = np.array(gas_load_info['gas_load_min'])
        self.gas_load_max = np.array(gas_load_info['gas_load_max'])
        self.gas_load_connection_index = np.array(gas_load_info['gas_load_connection_index'])
        # ----------------- gas generator info -----------------------
        self.gas_generator_num = gas_generator_info['gen_gas_num']
        self.gas_gen_index_power = np.array(gas_generator_info['gen_gas_index_power'])
        self.gas_gen_index_gas = np.array(gas_generator_info['gen_gas_index_gas'])
        self.gas_gen_capacity = gas_generator_info['gen_gas_capacity']
        # ----------------- gas well info ---------------------------
        self.gas_well_num = gas_well_info['gas_well_num']
        self.gas_well_connection_index = np.array(gas_well_info['gas_well_connection_index'])
        self.well_cost = gas_well_info['well_cost']
        # ------------------- gas transient START ---------------------------
        self.time_count_per_second = pipelineConfig['tc']  # tao / delta_h <= 1 / c
        self.pipeline_length = pipelineConfig['L']
        self.w_average = pipelineConfig['w_av']
        self.c1 = pipelineConfig['c1']
        self.demeter = pipelineConfig['dij']  # this is pipeline Demeter
        self.lamb_da = pipelineConfig['Lambda']  # this is fraction coefficient

        self.time_counts = pipelineConfig['Time_all']  # 暂态时间slot
        self.T_long = pipelineConfig['T_long']

        self.T = self.time_counts  # T is the alias of time all
        self.time_per_T = int(self.T / self.T_long)
        self.form_matrix()
        # config the model
        self.model = gurobi.Model()

    def generate_trans_coeff(self):
        tao = 1 / self.time_count_per_second  # this is delta t
        delta_h = self.c1 / self.time_count_per_second  # this is delta x
        A = self.demeter * self.demeter / 4 * 3.1415926
        c2 = self.c1 * self.c1
        tao2 = tao * tao
        h2 = delta_h * delta_h
        C = self.lamb_da * self.w_average / 2 / self.demeter

        self.node_counts = int(self.pipeline_length / delta_h)  # this is space node number

        self.C1 = 1 / tao2 + C / tao
        self.C2 = 2 / tao2 - 2 * c2 / h2 + C / tao
        self.C3 = -1 / tao2
        self.C4 = c2 / h2
        self.C5 = c2 / h2

        self.K_mp = A / c2 * delta_h / tao

        self.KP1 = 1 / (tao2 * c2) + self.lamb_da * self.w_average / (2 * self.demeter * c2 * tao)
        self.KP2 = -2 / h2 + 2 / (c2 * tao2) + C / (c2 * tao)
        self.KP3 = -1 / (c2 * tao2)
        self.KP4 = 1 / h2
        self.KP5 = 1 / h2

        self.KP_M_t_1 = delta_h / A / tao
        self.KP_M_t = -1 * (delta_h / A / tao + self.lamb_da * self.w_average * delta_h / 2 / self.demeter / A)
        self.K_p_ini_x = self.lamb_da * self.w_average * delta_h / 2 / self.demeter / A

        self.p_mt = np.zeros((self.time_counts, self.time_counts * 2), dtype=np.float)

        self.Matrix_ab = np.zeros((self.time_counts, self.node_counts, self.time_counts * 2), dtype=np.float)
        self.Matrix_cd = np.zeros((self.time_counts, self.node_counts, 2 * self.time_counts), dtype=np.float)

        self.A_ana = np.zeros((self.time_counts, self.time_counts), dtype=np.float)
        self.B_ana = np.zeros((self.time_counts, self.time_counts), dtype=np.float)
        self.C_ana = np.zeros((self.time_counts, self.time_counts), dtype=np.float)
        self.D_ana = np.zeros((self.time_counts, self.time_counts), dtype=np.float)

    def generate_matrix_ab(self):
        # construct for time of the first two layer:
        for i in range(self.node_counts):
            self.Matrix_ab[0, i] = [1] + [0] * (self.time_counts - 1) + [0] * self.time_counts
            self.Matrix_ab[1, i] = [0] + [1] + [0] * (self.time_counts - 2) + [0] * self.time_counts

        # construct for the node of the first node through all time layer:
        for i in range(self.time_counts):
            self.Matrix_ab[i, 0] = [0] * i + [1] + [0] * (self.time_counts - 1 - i) + [0] * self.time_counts

        for t in range(self.time_counts):
            if t == 0:
                self.p_mt[t] = [0] * (2 * self.time_counts)
                self.p_mt[t] = [0] * self.time_counts + [-1] + [1] + [0] * (self.time_counts - 2)
                #                               two time diff
                continue
            self.p_mt[t] = [0] * self.time_counts + [0] * (t - 1) + [-1] * 1 + [1] * 1 + [0] * (
                        self.time_counts - 1 - t)
            #      all zero for first time_counts   zero for pre t     two time diff         zero for post t

        # for time start from the three layer:
        for t in range(self.time_counts - 2):
            # for the node start from the second node:
            for n in range(self.node_counts - 2):
                self.Matrix_ab[t + 2][n + 1] = \
                    (self.C2 * self.Matrix_ab[t + 1][n + 1] +
                     self.C3 * self.Matrix_ab[t][n + 1] +
                     self.C4 * self.Matrix_ab[t + 1, n + 2] +
                     self.C5 * self.Matrix_ab[t + 1][n]) / self.C1
            # for the node of the last node:
            self.Matrix_ab[t + 2][self.node_counts - 1] = \
                self.Matrix_ab[t + 2][self.node_counts - 2] + self.K_mp * self.p_mt[t + 2]

        self.Matrix_ab[:, -1, self.time_counts: 2 * self.time_counts] = \
            self.Matrix_ab[:, -1, self.time_counts:2 * self.time_counts] * 4e3

        for t in range(self.time_counts):
            self.A_ana[t, :] = self.Matrix_ab[t][-1][0:self.time_counts]
            self.B_ana[t, :] = self.Matrix_ab[t][-1][self.time_counts:2 * self.time_counts]

    def generate_matrix_cd(self):
        for node in range(self.node_counts):
            self.Matrix_cd[0, node] = [- node * self.K_p_ini_x] + [0] * (self.time_counts - 1) + [1] + [0] * (
                        self.time_counts - 1)
            self.Matrix_cd[1, node] = [0] + [- node * self.K_p_ini_x] + \
                                      [0] * (self.time_counts - 2) + \
                                      [0] + [1] + \
                                      [0] * (self.time_counts - 2)  # t=1

        for t in range(self.time_counts):
            self.Matrix_cd[t, 0] = [0] * self.time_counts + [0] * t + [1] + [0] * (self.time_counts - 1 - t)  # n=0
        for t in range(self.time_counts - 2):  # fix t
            for n in range(self.node_counts - 2):
                self.Matrix_cd[t + 2, n + 1] = \
                    (self.KP2 * self.Matrix_cd[t + 1, n + 1] +
                     self.KP3 * self.Matrix_cd[t, n + 1] +
                     self.KP4 * self.Matrix_cd[t + 1, n + 2] +
                     self.KP5 * self.Matrix_cd[t + 1, n]) / self.KP1
            self.Matrix_cd[t + 2, self.node_counts - 1] = \
                self.Matrix_cd[t + 2, self.node_counts - 2] + \
                np.array([0] * (t + 2) + [1] + [0] * (2 * self.time_counts - t - 3)) * self.KP_M_t + \
                np.array(([0] * (t + 1) + [1] + [0] * (2 * self.time_counts - t - 2))) * self.KP_M_t_1

        self.Matrix_cd[:, -1, 0: self.time_counts] = self.Matrix_cd[:, -1, 0: self.time_counts] / 4e3

        for t in range(self.time_counts):
            self.C_ana[t, :] = self.Matrix_cd[t][-1][0:self.time_counts]
            self.D_ana[t, :] = self.Matrix_cd[t][-1][self.time_counts:2 * self.time_counts]

    def form_matrix(self):
        import os.path as path
        import pickle

        if not path.exists('abcd.pkl'):
            self.generate_trans_coeff()
            self.generate_matrix_ab()
            self.generate_matrix_cd()
            with open('abcd.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
                pickle.dump([self.Matrix_ab, self.Matrix_cd, self.node_counts], f)
        else:
            print('load matrix from file')
            with open('abcd.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
                self.Matrix_ab, self.Matrix_cd, self.node_counts = pickle.load(f)
                self.A_ana = self.Matrix_ab[:, -1, 0:self.time_counts]
                self.B_ana = self.Matrix_ab[:, -1, self.time_counts: 2 * self.time_counts]
                self.C_ana = self.Matrix_cd[:, -1, 0:self.time_counts]
                self.D_ana = self.Matrix_cd[:, -1, self.time_counts: 2 * self.time_counts]

    def build_original_model_var(self):
        # build power system variables
        def getTranVars(varLong):
            rows = 0
            if isinstance(varLong, gurobi.tupledict):
                rows = varLong.keys()[-1][0] + 1  # +1 is actual length
            elif isinstance(varLong, np.ndarray):
                rows = varLong.shape[0]
            trans = np.empty((rows, self.T), dtype=np.object)
            for row in range(rows):
                for t in range(self.T):
                    trans[row, t] = varLong[row, int(t / self.time_per_T)]
            return trans

        self.power_gen = tonp(self.model.addVars(
            self.gen_num, self.T_long,
            lb=[[self.gen_power_min[i]] * self.T_long for i in range(self.gen_num)],
            ub=[[self.gen_power_max[i]] * self.T_long for i in range(self.gen_num)],
            name='power_gene'))
        self.react_gen = tonp(self.model.addVars(
            self.gen_num, self.T_long,
            lb=[[self.gen_react_min[i]] * self.T_long for i in range(self.gen_num)],
            ub=[[self.gen_react_max[i]] * self.T_long for i in range(self.gen_num)],
            name='reactive_gene'))
        self.power_load = tonp(self.model.addVars(
            self.load_num, self.T_long,
            lb=self.load_power_min,  # [[self.load_power_min[i]] * self.T for i in range(self.load_num)],
            ub=self.load_power_max,  # [[self.load_power_max[i]] * self.T for i in range(self.load_num)],
            name='power_load'))
        self.react_load = tonp(self.model.addVars(
            self.load_num, self.T_long,
            lb=self.load_react_min,  # [[self.load_react_min[i]] * self.T for i in range(self.load_num)],
            ub=self.load_react_max,  # [[self.load_react_max[i]] * self.T for i in range(self.load_num)],
            name='react_load'))

        self.power_gen_trans = getTranVars(self.power_gen)
        self.react_gen_trans = getTranVars(self.react_gen)
        self.power_load_trans = getTranVars(self.power_load)
        self.react_load_trans = getTranVars(self.react_load)

        self.voltage_square = tonp(self.model.addVars(
            self.bus_num, self.T_long,
            lb=self.bus_voltage_min[0] * self.bus_voltage_min[0],
            ub=self.bus_voltage_max[0] * self.bus_voltage_max[0],
            name='bus_voltage_square'))
        self.line_current_square = tonp(self.model.addVars(
            self.line_num, self.T_long,
            ub=self.line_current_capacity[0] * self.line_current_capacity[0],
            name='line_current_square'))
        self.line_power_flow = tonp(self.model.addVars(
            self.line_num, self.T_long,
            lb=-1 * self.line_power_flow_max, ub=self.line_power_flow_max,
            name='line_power_flow'))
        self.line_react_flow = tonp(self.model.addVars(
            self.line_num, self.T_long,
            lb=-1 * self.line_react_flow_max, ub=self.line_react_flow_max,
            name='line_react_flow'))

        # build gas system variables
        self.start_pipeline_flow_trans = tonp(self.model.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_flow_min, ub=self.gas_flow_max,
            name='pipeline_gas_flow_of_start_port'
        ))
        self.end_pipeline_flow_trans = tonp(self.model.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_flow_min, ub=self.gas_flow_max,
            name='pipeline_gas_flow_of_end_port'
        ))
        self.start_pipeline_pressure_trans = tonp(self.model.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
            name='pipeline_gas_pressure_of_start_port'
        ))
        self.end_pipeline_pressure_trans = tonp(self.model.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
            name='pipeline_gas_pressure_of_end_port'
        ))

        ########
        self.node_pressure_trans = tonp(self.model.addVars(
            self.gas_node_num, self.T,
            lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
            name='gas_node_pressure_trans'
        ))
        ########

        # output of gas generator is T-long, but the reserve is real-time
        self.gas_generator = tonp(self.model.addVars(
            self.gas_generator_num, self.T_long,
            name='gas_generator_transient_output'
        ))
        self.gas_generator_tran_reserve_trans = tonp(self.model.addVars(
            self.gas_generator_num, self.T,
            name='gas_generator_transient_reserve'
        ))
        # gas load is T-long
        self.gas_load = self.model.addVars(
            self.gas_load_num, self.T_long,
            lb=self.gas_load_min,
            ub=self.gas_load_max,
            name='gas_load_transient_output'
        )
        # gas well gas output in real-time
        self.well_in_trans = tonp(self.model.addVars(
            self.gas_well_num, self.T,
            name='gas_well_transient_output'
        ))
        self.well_pressure = self.model.addVars(
            self.gas_well_num, self.T_long,
            name='gas_well_set_pressure'
        )
        self.gas_generator_trans = getTranVars(self.gas_generator)
        self.gas_load_trans = getTranVars(self.gas_load)
        self.well_pressure_trans = getTranVars(self.well_pressure)

    def buildBasePrimaryProblem(self):
        # ----------- 1. add Vars ---------------------------
        self.build_original_model_var()
        # ----------- 2.1 node power balance -----------------
        self.model.update()
        #
        print('start')
        for node in range(self.bus_num):
            for time in range(self.T_long):
                self.model.addConstr(
                    sum(self.power_gen[np.where(np.array(self.gen_conn_power_index) == node), time].flatten()) +
                    sum(self.line_power_flow[np.where(np.array(self.line_end_point) == node), time].flatten()) -
                    sum((self.line_resistance[np.where(np.array(self.line_end_point) == node)] *
                         self.line_current_square[np.where(np.array(self.line_end_point) == node), time]).flatten()) +
                    sum(self.gas_generator[np.where(np.array(self.gas_gen_index_power) == node), time].flatten())
                    ==
                    sum(self.power_load[np.where(np.array(self.load_index) == node), time].flatten()) +
                    sum(self.line_power_flow[np.where(np.array(self.line_start_point) == node), time].flatten()),
                    name='power_balance')
                self.model.addConstr(
                    sum(self.react_gen[np.where(np.array(self.gen_conn_power_index) == node), time].flatten()) +
                    sum(self.line_react_flow[np.where(np.array(self.line_end_point) == node), time].flatten()) -
                    sum((self.line_reactance[np.where(np.array(self.line_end_point) == node)] *
                         self.line_current_square[np.where(np.array(self.line_end_point) == node), time]).flatten())
                    ==
                    sum(self.react_load[np.where(np.array(self.load_index) == node), time].flatten()) +
                    sum(self.line_react_flow[np.where(np.array(self.line_start_point) == node), time].flatten()),
                    name='react_balance')



        # ----------- 2.2 line voltage drop ------------------
        for i in range(self.line_num):
            start_point, end_point = self.line_start_point[i], self.line_end_point[i]
            resistance, reactance = self.line_resistance[i], self.line_reactance[i]
            impedance_square = reactance * reactance + resistance * resistance
            for time in range(self.T_long):
                self.model.addConstr(
                    self.voltage_square[end_point, time] -
                    self.voltage_square[start_point, time]
                    ==
                    impedance_square * self.line_current_square[i, time] -
                    2 * (resistance * self.line_power_flow[i, time] +
                         reactance * self.line_react_flow[i, time]),
                    name='voltage_drop')
                self.model.addConstr(
                    self.line_power_flow[i, time] * self.line_power_flow[i, time] +
                    self.line_react_flow[i, time] * self.line_react_flow[i, time]
                    <=
                    self.line_current_square[i, time] * self.voltage_square[start_point, time],
                    name='flow_relax')

        # ====> Two port network of gas pipeline [gas system transient]
        # self.Matrix_ab[self.Matrix_ab <= 1e-8] = 0
        # self.Matrix_cd[self.Matrix_cd <= 1e-8] = 0
        def myDot2(matrix, syms, row, y_ranges):
            value = matrix[row][self.node_counts - 1][y_ranges]
            sym = syms[y_ranges]
            return sum(value * sym)

        # ====> 2.3 two port pipeline [gas system transient]
        # print('======> finish 7.1')
        import time
        t333 = time.time()
        x, y = np.nonzero(self.Matrix_ab[:, self.node_counts - 1, :])
        x2, y2 = np.nonzero(self.Matrix_cd[:, self.node_counts - 1, :])
        for t in range(self.T):
            y_range = y[np.where(x == t)]  # 0 1
            y_range2 = y2[np.where(x2 == t)]
            for line in range(self.gas_line_num):
                # 1              2
                # o ------------ o
                # start          end
                start_node = self.gas_pipeline_start_node[line]
                end_node = self.gas_pipeline_end_node[line]

                start_pipeline_flow = self.start_pipeline_flow_trans[line]
                end_pipeline_flow = self.end_pipeline_flow_trans[line]
                start_pipeline_pressure = self.node_pressure_trans[start_node]
                end_pipeline_pressure = self.node_pressure_trans[end_node]

                flow2pressure1 = np.hstack((end_pipeline_flow, start_pipeline_pressure))

                # print('=======> finish 7.1.1')
                self.model.addConstr(start_pipeline_flow[t] ==
                                     myDot2(self.Matrix_ab, flow2pressure1, t, y_range), name='c18')
                self.model.addConstr(end_pipeline_pressure[t] ==
                                     myDot2(self.Matrix_cd, flow2pressure1, t, y_range2), name='c17')

                # self.model.addConstr(start_pipeline_flow[t] ==
                #                      sum(self.Matrix_ab[t][self.node_counts - 1] * flow2pressure1), name='c18')
                # self.model.addConstr(end_pipeline_pressure[t] ==
                #                      sum(self.Matrix_cd[t][self.node_counts - 1] * flow2pressure1), name='c17')

        t2222 = time.time()
        print(t2222 - t333)
        # ====> 2.4 node gas balance [gas system transient]
        # print('======> finish 7.2')
        for t in range(self.T):
            for node in range(self.gas_node_num):
                self.model.addConstr(
                    sum(self.start_pipeline_flow_trans[np.where(self.gas_pipeline_start_node == node),   t].flatten()) +
                    sum(self.gas_load_trans[           np.where(self.gas_load_connection_index == node), t].flatten()) +
                    sum(self.gas_generator_trans[      np.where(self.gas_gen_index_gas == node),         t].flatten())
                    ==
                    sum(self.end_pipeline_flow_trans[  np.where(self.gas_pipeline_end_node == node),     t].flatten()) +
                    sum(self.start_pipeline_flow_trans[np.where(self.gas_pipeline_start_node == node)
                                                       if node in self.gas_well_connection_index
                                                       else tuple(np.array([])),                         t].flatten()),
                    name='c19')

        # ====> 2.6 gas generator transient [gas system transient]
        # print('======> finish 7.3')
        for t in range(self.T):
            for gen in range(self.gas_generator_num):
                self.model.addConstr(
                    self.gas_generator_trans[gen, t] + self.gas_generator_tran_reserve_trans[gen, t] <=
                    self.gas_gen_capacity[gen], name='c22'
                )
                self.model.addConstr(
                    self.gas_generator_trans[gen, t] - self.gas_generator_tran_reserve_trans[gen, t] >=
                    0, name='c31'
                )

        # ====> 2.7 node pressure specify [gas system transient]
        # print('======> finish 7.5')
        for t in range(self.T):
            for well in range(self.gas_well_num):
                node = self.gas_well_connection_index[well]

                self.model.addConstr(
                    self.node_pressure_trans[node, t] ==
                    self.PRESSURE_CONSTANT, name='c26'
                )



        objs = []
        for load in range(self.load_num):
            for time in range(self.T_long):
                objs.append((self.power_load[load, time] - (self.load_power_max[load, time] + self.load_power_min[load, time]) / 2) *
                            (self.power_load[load, time] - (self.load_power_max[load, time] + self.load_power_min[load, time]) / 2))
                objs.append((self.react_load[load, time] - (self.load_react_max[load, time] + self.load_react_min[load, time]) / 2) *
                            (self.react_load[load, time] - (self.load_react_max[load, time] + self.load_react_min[load, time]) / 2))
        for gen in range(self.gen_num):
            for time in range(self.T_long):
                objs.append(self.gen_cost_a[gen] * self.power_gen[gen, time] * self.power_gen[gen, time])
                objs.append(self.gen_cost_b[gen] * self.power_gen[gen, time])
                objs.append(self.gen_cost_c[gen])
        for well in range(self.gas_well_num):
            for time in range(self.T):
                objs.append(self.start_pipeline_flow_trans[4, time] * self.well_cost / self.time_counts)
        for load in range(self.gas_load_num):
            for time in range(self.T_long):
                objs.append((self.gas_load[load, time] - (self.gas_load_min[load, time] + self.gas_load_max[load, time]) / 2) *
                            (self.gas_load[load, time] - (self.gas_load_min[load, time] + self.gas_load_max[load, time]) / 2))
        # --------------- 3 set objective ----------------
        self.model.setObjective(sum(objs))
        # self.model.setParam('BarHomogeneous', 1)



if __name__ == '__main__':
    pg = PowerGas(*getConfig())

    pg.buildBasePrimaryProblem()
    pg.model.optimize()