import matplotlib.pyplot as plt
import numpy as np
import gurobipy as gurobi
from config.powergas1 import getConfig
import time

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
            if vars_gur[i] is not None:
                re[i] = vars_gur[i].getAttr('X')
            else:
                re[i] = 0
        return re
    if dim == 2:
        re = np.empty(tuple(shape), dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                if vars_gur[i, j] is not None:
                    re[i, j] = vars_gur[i, j].getAttr('X')
                else:
                    re[i, j] = 0
        return re
    if dim == 3:
        re = np.empty(shape, dtype=np.float)
        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    if vars_gur[i, j, k] is not None:
                        re[i, j, k] = vars_gur[i, j, k].getAttr('X')
                    else:
                        re[i, j, k] = 0
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

def getVarsDiff(expression):
        assert (isinstance(expression, gurobi.LinExpr))

        expr = expression
        varSize = expr.size()
        allVars = []
        allCoeff = []

        for index in range(varSize):
            allVars.append(expr.getVar(index))
            allCoeff.append(expr.getCoeff(index))

        constant = expr.getConstant()

        var_value = {}
        for index, var in enumerate(allVars):
            if var not in var_value:
                var_value[var] = allCoeff[index]
            else:
                var_value[var] += allCoeff[index]

        allVars = list(var_value.keys())
        allCoeff = list(var_value.values())
        return allVars, allCoeff, constant
INF = gurobi.GRB.INFINITY

LOG_OUTPUT = True
TIME_LONG_TO_FRE = 2

# TODO:
#  1, original problem solve
#  2, dual problem solve
#  3, add constraints algorithm

class PowerGas:
    def __init__(self, powerConfig, gasConfig, pipelineConfig):
        self.picture_index = 1
        self.file_index = 0
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
        self.T_fre = self.T_long * TIME_LONG_TO_FRE

        self.T = self.time_counts  # T is the alias of time all
        self.time_per_T = int(self.T / self.T_long)
        self.form_matrix()
        # config the model
        self.model = gurobi.Model()
        assert (self.T % self.T_long) == 0
        assert (self.T % self.T_fre) == 0

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

        # if not path.exists('abcd.pkl'):
        #     self.generate_trans_coeff()
        #     self.generate_matrix_ab()
        #     self.generate_matrix_cd()
        #     with open('abcd.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
        #         pickle.dump([self.Matrix_ab, self.Matrix_cd, self.node_counts], f)
        # else:
        #     print('load matrix from file')
        #     with open('abcd.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        #         self.Matrix_ab, self.Matrix_cd, self.node_counts = pickle.load(f)
        #         self.A_ana = self.Matrix_ab[:, -1, 0:self.time_counts]
        #         self.B_ana = self.Matrix_ab[:, -1, self.time_counts: 2 * self.time_counts]
        #         self.C_ana = self.Matrix_cd[:, -1, 0:self.time_counts]
        #         self.D_ana = self.Matrix_cd[:, -1, self.time_counts: 2 * self.time_counts]

        self.generate_trans_coeff()
        self.generate_matrix_ab()
        self.generate_matrix_cd()
        with open('abcd.pkl', 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([self.Matrix_ab, self.Matrix_cd, self.node_counts], f)

        self.generate_network_matrix()
        self.testMatrixHaHa()
    def testMatrixHaHa(self):
        flow_load = np.vstack((np.ones((int(self.T / 3), 1)),
                               np.ones((int(self.T / 2) - int(self.T / 3), 1)) * 3,
                               np.ones((self.T - int(self.T / 2), 1)) * 1
                               ))
        gas_generator = np.ones((self.T, 1)) * 2
        pressure_well = np.ones((self.T, 1)) * 100

        len1 = flow_load.shape[0]
        len2 = gas_generator.shape[0]
        len3 = pressure_well.shape[0]

        MldPsr = np.vstack((
            flow_load.reshape(-1, 1),
            gas_generator.reshape(-1, 1),
            pressure_well.reshape(-1, 1)
        ))

        flow_source = np.zeros_like(flow_load)
        pressure_load = np.zeros_like(flow_load)
        pressure_generator = np.zeros_like(flow_load)
        MsrPld = np.vstack((
            flow_source        .reshape(-1, 1),
            pressure_load      .reshape(-1, 1),
            pressure_generator .reshape(-1, 1)
        ))

        YY = self.YY
        temp = YY.dot(MldPsr)

        flow_source = temp[:len1, 0] * -1
        pressure_load = temp[len1: len1 + len2, 0]
        pressure_generator = temp[len1+len2: , 0]

        import matplotlib.pyplot as plt
        plt.subplot(3, 1, 1)
        plt.plot(flow_source)

        plt.subplot(3, 1, 2)
        plt.plot(pressure_load)

        plt.subplot(3, 1, 3)
        plt.plot(pressure_generator)
        plt.show()
    def generate_network_matrix(self):
        # generate node-branch matrix aka Ain
        self.node_branch_end = np.zeros((self.gas_node_num * self.T, self.gas_line_num * self.T))
        for node in range(self.gas_node_num):
            index = np.where(np.array(self.gas_pipeline_end_node) == node)[0]
            row = np.zeros((self.T, self.T * self.gas_line_num))

            for line in range(self.gas_line_num):
                if line in index:
                    row[:, line*self.T:(line+1)*self.T] = np.eye(self.T, self.T)

            self.node_branch_end[node*self.T:(node+1)*self.T, :] = row

        # generate node-branch matrix aka Aout
        self.node_branch_start = np.zeros((self.gas_node_num * self.T, self.gas_line_num * self.T))
        for node in range(self.gas_node_num):
            index = np.where(np.array(self.gas_pipeline_start_node) == node)[0]
            row = np.zeros((self.T, self.T * self.gas_line_num))

            for line in range(self.gas_line_num):
                if line in index:
                    row[:, line*self.T:(line+1)*self.T] = np.eye(self.T, self.T)

            self.node_branch_start[node*self.T:(node+1)*self.T, :] = row

        # generate line-start-node pressure index matrix, aka Apn1
        self.branch_node_start = np.zeros((self.gas_line_num * self.T, self.gas_node_num * self.T))
        for line in range(self.gas_line_num):
            start_node = self.gas_pipeline_start_node[line]
            row = np.zeros((self.T, self.T * self.gas_node_num))
            row[:, start_node * self.T : (start_node + 1) * self.T] = np.eye(self.T, self.T)

            self.branch_node_start[line*self.T:(line+1)*self.T, :] = row

        # generate line-end-node pressure index matrix, aka Apn2
        self.branch_node_end = np.zeros((self.gas_line_num * self.T, self.gas_node_num * self.T))
        for line in range(self.gas_line_num):
            end_node = self.gas_pipeline_end_node[line]
            row = np.zeros((self.T, self.T * self.gas_node_num))
            row[:, end_node * self.T : (end_node + 1) * self.T] = np.eye(self.T, self.T)

            self.branch_node_end[line*self.T:(line+1)*self.T, :] = row

        A = self.A_ana
        B = self.B_ana
        C = self.C_ana
        Cinv = np.linalg.inv(C)
        D = self.D_ana
        Ain = self.node_branch_end
        Aout = self.node_branch_start
        Apn1 = self.branch_node_start
        Apn2 = self.branch_node_end

        # get matrix Y
        first = np.block([-1 * Aout, Ain])
        second = np.block([
            [B - A.dot(Cinv).dot(D), A.dot(Cinv)],
            [-1 * Cinv.dot(D),       Cinv]])

        import scipy.linalg as lin
        from numpy.linalg import inv as inv
        n = self.gas_line_num

        def getFull(matrix):
            mat = matrix.copy()
            for i in range(n - 1):
                mat = lin.block_diag(mat, matrix)
            return mat

        full_second = np.block([
            [getFull(B - A.dot(Cinv).dot(D)), getFull(A.dot(Cinv))],
            [getFull(-1 * Cinv.dot(D)),       getFull(Cinv)]
        ])

        third = np.block([
            [Apn1],
            [Apn2]
        ])
        Y = first.dot(full_second).dot(third)

        # so we assume the first one node is well node
        #              then is internal node
        #              last is load node
        internal_count = self.gas_node_num - self.gas_load_num - \
                         self.gas_generator_num - self.gas_well_num
        # load_count = 2

        Y11 = Y[0:self.T, 0:self.T]
        Y12 = Y[0:self.T, self.T: self.T + self.T * internal_count]
        Y13 = Y[0:self.T, self.T + self.T * internal_count : ]

        Y21 = Y[self.T:self.T + self.T * internal_count, 0:self.T]
        Y22 = Y[self.T:self.T + self.T * internal_count, self.T: self.T + self.T * internal_count]
        Y23 = Y[self.T:self.T + self.T * internal_count, self.T + self.T * internal_count : ]

        Y31 = Y[self.T + self.T * internal_count:, 0:self.T]
        Y32 = Y[self.T + self.T * internal_count:, self.T: self.T + self.T * internal_count]
        Y33 = Y[self.T + self.T * internal_count:, self.T + self.T * internal_count : ]

        Ymm = (Y13 - Y12.dot(inv(Y22)).dot(Y23)).dot(
               inv(Y33 - Y32.dot(inv(Y22)).dot(Y23)) )

        Ymp = (Y11 - Y12.dot(inv(Y22)).dot(Y21)) - \
              ((Y13 - Y12.dot(inv(Y22)).dot(Y23)).dot(
               inv(Y33 - Y32.dot(inv(Y22)).dot(Y23))).dot(
               Y31 - Y32.dot(inv(Y22)).dot(Y21)))

        Ypm = inv(Y33 - Y32.dot(inv(Y22)).dot(Y23))

        Ypp = -1 * (inv(Y33 - Y32.dot(inv(Y22)).dot(Y23))).dot(
               Y31 - Y32.dot(inv(Y22)).dot(Y21))

        self.Ymm = Ymm
        self.Ymp = Ymp
        self.Ypm = Ypm
        self.Ypp = Ypp

        self.YY = np.block([
            [self.Ymm, self.Ymp],
            [self.Ypm, self.Ypp]
        ])




    def build_original_model_var_network(self):
        # build power system variables
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

        # self.power_gen_trans = self.getTranVars(self.power_gen)
        # self.react_gen_trans = self.getTranVars(self.react_gen)
        # self.power_load_trans = self.getTranVars(self.power_load)
        # self.react_load_trans = self.getTranVars(self.react_load)

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



        self.gas_generator = tonp(self.model.addVars(
            self.gas_generator_num, self.T_long,
            name='gas_generator_transient_output'
        ))
        self.gas_generator_reserve_up = tonp(self.model.addVars(
            self.gas_generator_num, self.T_long,
            name='gas_generator_transient_reserve'
        ))
        self.gas_generator_reserve_down = tonp(self.model.addVars(
            self.gas_generator_num, self.T_long,
            name='gas_generator_transient_reserve'
        ))

        # gas load is T-long
        # 气负荷 小时级
        self.gas_load = self.model.addVars(
            self.gas_load_num, self.T_long,
            lb=self.gas_load_min,
            ub=self.gas_load_max,
            name='gas_load_transient_output'
        )
        # 气井出力
        # gas well gas output in real-time
        self.well_in_trans = tonp(self.model.addVars(
            self.gas_well_num, self.T,
            name='gas_well_transient_output'
        ))
        self.well_pressure = tonp(self.model.addVars(
            self.gas_well_num, self.T_long,
            name='gas_well_set_pressure'
        ))

        self.gas_flow_source = tonp(self.model.addVars(
            1, self.T,
            lb=-1*INF, ub=0,
            name='gas_flow_source'
        ))
        self.gas_pressure_load = tonp(self.model.addVars(
            self.gas_load_num, self.T
        ))
        self.gas_pressure_generator = tonp(self.model.addVars(
            self.gas_generator_num, self.T
        ))

        self.gas_generator_trans = self.getTranVars(self.gas_generator)
        self.gas_load_trans = self.getTranVars(self.gas_load)
        self.well_pressure_trans = self.getTranVars(self.well_pressure)
    def buildBasePrimaryProblemNetwork(self):
        # ----------- 1. add Vars ---------------------------
        self.build_original_model_var_network()
        # ----------- 2.1 node power balance -----------------
        self.model.update()

        print('start')
        for node in range(self.bus_num):
            for t in range(self.T_long):
                self.model.addConstr(
                    sum(self.power_gen[np.where(np.array(self.gen_conn_power_index) == node), t].flatten()) +
                    sum(self.line_power_flow[np.where(np.array(self.line_end_point) == node), t].flatten()) -
                    sum((self.line_resistance[np.where(np.array(self.line_end_point) == node)] *
                         self.line_current_square[np.where(np.array(self.line_end_point) == node), t]).flatten()) +
                    sum(self.gas_generator[np.where(np.array(self.gas_gen_index_power) == node), t].flatten())
                    ==
                    sum(self.power_load[np.where(np.array(self.load_index) == node), t].flatten()) +
                    sum(self.line_power_flow[np.where(np.array(self.line_start_point) == node), t].flatten()),
                    name='power_balance')
                self.model.addConstr(
                    sum(self.react_gen[np.where(np.array(self.gen_conn_power_index) == node), t].flatten()) +
                    sum(self.line_react_flow[np.where(np.array(self.line_end_point) == node), t].flatten()) -
                    sum((self.line_reactance[np.where(np.array(self.line_end_point) == node)] *
                         self.line_current_square[np.where(np.array(self.line_end_point) == node), t]).flatten())
                    ==
                    sum(self.react_load[np.where(np.array(self.load_index) == node), t].flatten()) +
                    sum(self.line_react_flow[np.where(np.array(self.line_start_point) == node), t].flatten()),
                    name='react_balance')

        # ----------- 2.2 line voltage drop ------------------
        for i in range(self.line_num):
            start_point, end_point = self.line_start_point[i], self.line_end_point[i]
            resistance, reactance = self.line_resistance[i], self.line_reactance[i]
            impedance_square = reactance * reactance + resistance * resistance
            for t in range(self.T_long):
                self.model.addConstr(
                    self.voltage_square[end_point, t] -
                    self.voltage_square[start_point, t]
                    ==
                    impedance_square * self.line_current_square[i, t] -
                    2 * (resistance * self.line_power_flow[i, t] +
                         reactance * self.line_react_flow[i, t]),
                    name='voltage_drop')
                self.model.addConstr(
                    self.line_power_flow[i, t] * self.line_power_flow[i, t] +
                    self.line_react_flow[i, t] * self.line_react_flow[i, t]
                    <=
                    self.line_current_square[i, t] * self.voltage_square[start_point, t],
                    name='flow_relax')

        # ====> Two port network of gas pipeline [gas system transient]







        # build network constraints
        MldPsr = np.vstack((
            self.gas_load_trans             .reshape(-1, 1),         #
            self.gas_generator_trans        .reshape(-1, 1),         # this is var
            self.well_pressure_trans        .reshape(-1, 1)          # this is known value
        ))

        MsrPld = np.vstack((
            self.gas_flow_source           .reshape(-1, 1),
            self.gas_pressure_load         .reshape(-1, 1),
            self.gas_pressure_generator    .reshape(-1, 1)       # all these are vars
        ))
        # byPassIndex = self.append_flow_source.reshape(-1, 1).shape[0]

        if LOG_OUTPUT:
            print('----add prim constr Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count):
            if row % 100 == 0:
                if LOG_OUTPUT:
                    print('------- add prim constraint' + str(row) + '----------')
            self.model.addConstr(
                MsrPld[row, 0] ==
                (self.YY[row, :]).dot(MldPsr)[0]
            )


        if LOG_OUTPUT:
            print('----add prim constr robust Step 3----------------')
        # build pressure limit
        node_collection = self.gas_pressure_load.flatten().tolist() + \
                          self.gas_pressure_generator.flatten().tolist()
        for load in node_collection:
            self.model.addConstr(
                load <= self.gas_node_pressure_max
            )
            self.model.addConstr(
                load >= self.gas_node_pressure_min
            )

        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                self.model.addConstr(
                    self.gas_generator[gen, t] + self.gas_generator_reserve_up[gen, t] <= self.gas_gen_capacity[gen]
                )
                self.model.addConstr(
                    self.gas_generator[gen, t] - self.gas_generator_reserve_down[gen, t] >= 0
                )

        for t in range(self.T):
            self.model.addConstr(
                self.well_pressure_trans[0, t] == self.PRESSURE_CONSTANT
            )

        objs = []

        for load in range(self.load_num):
            for t in range(self.T_long):
                objs.append((self.power_load[load, t] - (self.load_power_max[load, t] + self.load_power_min[load, t]) / 2) *
                            (self.power_load[load, t] - (self.load_power_max[load, t] + self.load_power_min[load, t]) / 2))
                objs.append((self.react_load[load, t] - (self.load_react_max[load, t] + self.load_react_min[load, t]) / 2) *
                            (self.react_load[load, t] - (self.load_react_max[load, t] + self.load_react_min[load, t]) / 2))

        for gen in range(self.gen_num):
            for t in range(self.T_long):
                objs.append(self.gen_cost_a[gen] * self.power_gen[gen, t] * self.power_gen[gen, t])
                objs.append(self.gen_cost_b[gen] * self.power_gen[gen, t])
                objs.append(self.gen_cost_c[gen])

        for well in range(self.gas_well_num):
            for t in range(self.T):
                objs.append(-1 * self.gas_flow_source[well, t] * self.well_cost / self.time_counts)

        for load in range(self.gas_load_num):
            for t in range(self.T_long):
                objs.append((self.gas_load[load, t] - (self.gas_load_min[load, t] + self.gas_load_max[load, t]) / 2) *
                            (self.gas_load[load, t] - (self.gas_load_min[load, t] + self.gas_load_max[load, t]) / 2))

        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                objs.append(-1 * self.gas_generator_reserve_up[gen, t] - self.gas_generator_reserve_down[gen, t])
        # --------------- 3 set objective ----------------
        # self.model.setParam('OutputFlag', 0)
        # self.model.setParam('BarHomogeneous', 1)
        # self.model.setParam('Method', 1)
        self.model.setObjective(sum(objs))









    # [USE] this is actually use network-form feasible problem
    def getTranVars(self, varLong, time_interval=0):
        if time_interval == 0:
            time_interval = self.time_per_T
        rows = 0
        if isinstance(varLong, gurobi.tupledict):
            rows = varLong.keys()[-1][0] + 1  # +1 is actual length
        elif isinstance(varLong, np.ndarray):
            rows = varLong.shape[0]
        trans = np.empty((rows, self.T), dtype=np.object)
        for row in range(rows):
            for t in range(self.T):
                trans[row, t] = varLong[row, int(t/time_interval)]
        return trans
    def addRelaxConstraints(self, v, z, dual, BIGM, model):
        model.addConstr(-1 * BIGM * z <= v)
        model.addConstr(v <= 0)
        model.addConstr(-1 * BIGM * (1 - z) <= (dual - v))
        model.addConstr((dual - v) <= 0)
    def buildRobustVarsNetwork(self, model):
        # add power system vars
        self.robust_power_gen_base = to_value(self.power_gen)      # this is pre known
        self.robust_react_gen_base = to_value(self.react_gen)      # this is known
        self.robust_power_load_base = to_value(self.power_load)    # this is known
        self.robust_react_load_base = to_value(self.react_load)    # this is known

        self.robust_voltage_square = tonp(model.addVars(self.bus_num, self.T_fre))
        self.robust_line_current_square = tonp(model.addVars(self.line_num, self.T_fre))
        self.robust_line_power_flow = tonp(model.addVars(self.line_num, self.T_fre))
        self.robust_line_react_flow = tonp(model.addVars(self.line_num, self.T_fre))

        self.robust_aux1 = tonp(model.addVars(self.line_num, self.T_fre, lb=-1*INF, ub=INF))
        self.robust_aux2 = tonp(model.addVars(self.line_num, self.T_fre, lb=-1*INF, ub=INF))

        # add gas system vars
        self.robust_pressure_well              = to_value(self.well_pressure_trans[0])                # this is known source pressure

        self.robust_gas_generator_base         = to_value(self.gas_generator)                         # this is known generator load
        self.robust_gas_generator_reserve_up   = to_value(self.gas_generator_reserve_up)              # this is known generator reserve
        self.robust_gas_generator_reserve_down = to_value(self.gas_generator_reserve_down)            # this is known generator reserve

        self.robust_gas_generator              = tonp(model.addVars(self.gas_generator_num, self.T_fre))
        self.robust_gas_generator_trans        = self.getTranVars(self.robust_gas_generator, time_interval=int(self.T/self.T_fre))

        self.robust_flow_load                  = to_value(self.gas_load_trans)                                        # this is known load gas flow

        self.robust_flow_source                = tonp(model.addVars(self.gas_well_num, self.T))                          # this is unknown source gas flow

        self.robust_pressure_load              = tonp(model.addVars(self.gas_load_num, self.T))                        # this is unknown load/gen node pressure
        self.robust_pressure_generator         = tonp(model.addVars(self.gas_generator_num, self.T))                        # this is unknown load/gen node pressure

        # self.zUp                               = tonp(model.addVars(self.gas_generator_num, self.T_fre, vtype=gurobi.GRB.BINARY))
        # self.zDown                             = tonp(model.addVars(self.gas_generator_num, self.T_fre, vtype=gurobi.GRB.BINARY))

    def feasibleProblemNetworkPreset(self):
        # build robust var
        self.robustModel = gurobi.Model()
        self.dualGen = GenDual2(self.robustModel)

        rM = self.robustModel
        self.buildRobustVarsNetwork(rM)
    def feasibleProblemNetworkBuild(self, robustModel):
        dG = self.dualGen
        sCollection = []

        # build power system constraints
        for node in range(self.bus_num):
            for t in range(self.T_fre):
                gen_index = np.where(np.array(self.gen_conn_power_index) == node)
                line_end_index = np.where(np.array(self.line_end_point) == node)
                line_start_index = np.where(np.array(self.line_start_point) == node)
                gas_index = np.where(np.array(self.gas_gen_index_power) == node)
                power_load_index = np.where(np.array(self.load_index) == node)

                dG.addConstr(
                    sum(self.robust_power_gen_base[gen_index, int(t/TIME_LONG_TO_FRE)].flatten()) +
                    sum(self.robust_line_power_flow[line_end_index, t].flatten()) -
                    sum((self.line_resistance[line_end_index] *
                         self.robust_line_current_square[line_end_index, t]).flatten()) +
                    sum(self.robust_gas_generator[gas_index, t].flatten()),

                    sum(self.robust_power_load_base[power_load_index, int(t/TIME_LONG_TO_FRE)].flatten()) +
                    sum(self.robust_line_power_flow[line_start_index, t].flatten()),
                    sense='==',
                    name='robust_power_balance')
                dG.addConstr(
                    sum(self.robust_react_gen_base[gen_index, int(t/TIME_LONG_TO_FRE)].flatten()) +
                    sum(self.robust_line_react_flow[line_end_index, t].flatten()) -
                    sum((self.line_reactance[line_end_index] *
                         self.robust_line_current_square[line_end_index, t]).flatten()),

                    sum(self.robust_react_load_base[power_load_index, int(t/TIME_LONG_TO_FRE)].flatten()) +
                    sum(self.robust_line_react_flow[line_start_index, t].flatten()),
                    sense='==',
                    name='robust_react_balance')

        # ----------- 2.2 line voltage drop ------------------
        for i in range(self.line_num):
            start_point, end_point = self.line_start_point[i], self.line_end_point[i]
            resistance, reactance = self.line_resistance[i], self.line_reactance[i]
            impedance_square = reactance * reactance + resistance * resistance
            for t in range(self.T_fre):
                dG.addConstr(
                    self.robust_voltage_square[end_point, t] -
                    self.robust_voltage_square[start_point, t],
                    impedance_square * self.robust_line_current_square[i, t] -
                    2 * (resistance * self.robust_line_power_flow[i, t] +
                         reactance * self.robust_line_react_flow[i, t]),
                    sense='==',
                    name='voltage_drop')
                dG.addConstr(
                    self.robust_aux1[i, t],
                    (self.robust_line_current_square[i, t] -
                     self.robust_voltage_square[start_point, t]) / 2,
                    sense='=='
                )
                dG.addConstr(
                    self.robust_aux2[i, t],
                    (self.robust_line_current_square[i, t] +
                     self.robust_voltage_square[start_point, t]) / 2,
                    sense='=='
                )
                dG.addConstrQ(
                    self.robust_line_power_flow[i, t] * self.robust_line_power_flow[i, t],
                    self.robust_line_react_flow[i, t] * self.robust_line_react_flow[i, t],
                    self.robust_aux1[i, t] * self.robust_aux1[i, t],
                    self.robust_aux2[i, t] * self.robust_aux2[i, t],
                    sense='<='
                )

        for node in range(self.bus_num):
            for t in range(self.T_fre):
                sMin = robustModel.addVar()
                sPlus = robustModel.addVar()
                sCollection.append(sMin)
                sCollection.append(sPlus)

                dG.addConstr(self.robust_voltage_square[node, t] - sMin,
                             self.bus_voltage_max[0] * self.bus_voltage_max[0],
                             '<=')
                dG.addConstr(self.robust_voltage_square[node, t] + sPlus,
                             self.bus_voltage_min[0] * self.bus_voltage_min[0],
                             '>=')
                dG.addConstr(
                    sMin, 0, sense='>='
                )
                dG.addConstr(
                    sPlus, 0, sense='>='
                )
                dG.addConstr(
                    self.robust_voltage_square[node, t], 0, sense='>='
                )

        for line in range(self.line_num):
            for t in range(self.T_fre):
                dG.addConstr(self.robust_line_current_square[line, t],
                             0,
                             '>=')
                dG.addConstr(self.robust_line_current_square[line, t],
                             self.line_current_capacity[0] * self.line_current_capacity[0],
                             '<=')

                dG.addConstr(self.robust_line_power_flow[line, t],
                             -1 * self.line_power_flow_max,
                             '>=')
                dG.addConstr(self.robust_line_power_flow[line, t],
                             self.line_power_flow_max,
                             '<=')

                dG.addConstr(self.robust_line_react_flow[line, t],
                             -1 * self.line_react_flow_max,
                             '>=')
                dG.addConstr(self.robust_line_react_flow[line, t],
                             self.line_react_flow_max,
                             '<=')


        # build network constraints
        MldPsr = np.vstack((
            self.robust_flow_load           .reshape(-1, 1),         # this is known value
            self.robust_gas_generator_trans       .reshape(-1, 1),         # this is var
            self.robust_pressure_well       .reshape(-1, 1)          # this is known value
        ))
        start_index = self.robust_flow_load.reshape(-1, 1).shape[0]
        middle_index = self.robust_flow_load.reshape(-1, 1).shape[0] + \
                       self.robust_gas_generator_trans.reshape(-1, 1).shape[0]

        MsrPld = np.vstack((
            self.robust_flow_source        .reshape(-1, 1),
            self.robust_pressure_load      .reshape(-1, 1),
            self.robust_pressure_generator .reshape(-1, 1)       # all these are vars
        ))

        byPassIndex = self.robust_flow_source.reshape(-1, 1).shape[0]

        # because of the O^2 complex of index search , so we pre-set variable index
        dG.preset_var_index((
            self.robust_gas_generator.reshape(-1, 1).flatten().tolist(),
            # self.robust_flow_source.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist()
        ))
        MAGIC_LIST = np.arange(self.robust_gas_generator.reshape(-1, 1).shape[0]).tolist()

        map_value_index = 0
        index_m = {}
        robustModel.update()
        for vars_preset in (
            self.robust_gas_generator.reshape(-1, 1).flatten().tolist(),
            # self.robust_flow_source.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist()
        ):
            vars_list = vars_preset
            for var in vars_list:
                index_m[var] = map_value_index
                map_value_index = map_value_index + 1

        if LOG_OUTPUT:
            print('----generate robust Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count - byPassIndex):
            row = row + byPassIndex
            if row % 100 == 0:
                if LOG_OUTPUT:
                    print('----generate robust Step 1 ' + str(row) + ' ----------------')
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                sense='=='
            )
            continue

            # coeffInput = [ -1]
            # coeff = self.YY[row, start_index:middle_index]
            # coeff = coeff.reshape((self.T_fre, -1))
            # coeffInput.extend( (1 * np.sum(coeff, axis=1)).tolist())
            #
            # varInput =   [ MsrPld[row, 0] ]
            # varInput.extend( self.robust_gas_generator.flatten().tolist() )
            #
            # varIndex =   [ index_m[MsrPld[row, 0]]]
            # varIndex. extend( MAGIC_LIST )
            #
            # constantInput = 1 * (
            #     self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
            #     self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            # )
            # assert(len(coeffInput) == len(varIndex))
            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     '<=',
            #     name='xx',
            #     appendIndex=False,
            #     coeffInput=coeffInput,
            #     varInput=varInput,
            #     constantInput=constantInput,
            #     varIndex=varIndex,
            #     showTime=False
            # )
            #
            # coeffInput = [ 1]
            # coeff = self.YY[row, start_index:middle_index]
            # coeff = coeff.reshape((self.T_fre, -1))
            # coeffInput.extend( (-1 * np.sum(coeff, axis=1)).tolist())
            # varInput =   [ MsrPld[row, 0] ]
            # varInput.extend( self.robust_gas_generator.flatten().tolist() )
            #
            # varIndex =   [ index_m[MsrPld[row, 0]]]
            # varIndex.extend( MAGIC_LIST )
            #
            # constantInput = -1 * (
            #     self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
            #     self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            # )
            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     '>=',
            #     name='xx',
            #     appendIndex=False,
            #     coeffInput=coeffInput,
            #     varInput=varInput,
            #     constantInput=constantInput,
            #     varIndex=varIndex,
            #     showTime=False
            # )

        if LOG_OUTPUT:
            print('----generate robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_fre):
                dG.addConstr(
                    self.robust_gas_generator[gen, t],
                    0,
                    sense='==',
                    appendIndex=True
                )

        if LOG_OUTPUT:
            print('----generate robust Step 3----------------')
        # build pressure limit

        node_collection = self.robust_pressure_load.flatten().tolist() + \
                          self.robust_pressure_generator.flatten().tolist()
        for load in node_collection:
            sMin = robustModel.addVar()
            sPlus = robustModel.addVar()
            sCollection.append(sMin)
            sCollection.append(sPlus)

            dG.addConstr(
                load - sMin, self.gas_node_pressure_max, sense='<='
            )
            dG.addConstr(
                load + sPlus, self.gas_node_pressure_min, sense='>='
            )

            dG.addConstr(
                load, 0, sense='>='
            )
            dG.addConstr(
                sMin, 0, sense='>='
            )
            dG.addConstr(
                sPlus, 0, sense='>='
            )

        if LOG_OUTPUT:
            print('----generate robust Step 3.1 add objective----------------')
        # build primitive objective
        dG.addObjectiveMin(gurobi.quicksum(sCollection))

        if LOG_OUTPUT:
            print('----generate robust Step 3.2 get dual problem----------------')
        # add additional objective
        dualModel = dG.getDual()

        if LOG_OUTPUT:
            print('----generate robust Step 3.3 add additional objective----------------')
        dualIndexList = dG.robustIndex
        dualIndexCount = 0
        BIGM = 1e7
        oobbjj = []


        self.zUp = tonp(dualModel.addVars(self.gas_generator_num, self.T_fre, vtype=gurobi.GRB.BINARY))
        self.zDown = tonp(dualModel.addVars(self.gas_generator_num, self.T_fre, vtype=gurobi.GRB.BINARY))

        self.zGreatV1 = np.empty((self.gas_generator_num, self.T_fre), dtype=np.object)
        self.zGreatV2 = np.empty((self.gas_generator_num, self.T_fre), dtype=np.object)
        self.zLessV1 = np.empty((self.gas_generator_num, self.T_fre), dtype=np.object)
        self.zLessV2 = np.empty((self.gas_generator_num, self.T_fre), dtype=np.object)


        for gen in range(self.gas_generator_num):
            for t in range(self.T_fre):    # per equal constraints
                dual1 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1
                dual2 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1

                assert not dual1.sameAs(dual2)

                zUp = self.zUp[gen, t]
                zDown = self.zDown[gen, t]
                dualModel.addConstr(zUp + zDown <= 1)

                # for great part
                gasBase = self.robust_gas_generator_base[gen, int(t/TIME_LONG_TO_FRE)]
                gasUpper = self.robust_gas_generator_reserve_up[gen, int(t/TIME_LONG_TO_FRE)]
                gasDown = self.robust_gas_generator_reserve_down[gen, int(t/TIME_LONG_TO_FRE)]
                # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
                v1 = dualModel.addVar(lb=-1*INF, ub=0)
                v2 = dualModel.addVar(lb=-1*INF, ub=0)
                self.zGreatV1[gen, t] = v1
                self.zGreatV2[gen, t] = v2
                # linear     v1 = zUp * dual1
                self.addRelaxConstraints(v1, zUp, dual1, BIGM, dualModel)
                # linear    v2 = zDown * dual1
                self.addRelaxConstraints(v2, zDown, dual1, BIGM, dualModel)
                oobbjj.append(dual1 * gasBase + gasUpper * v1 - gasDown * v2)

                # for less part
                # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
                v1 = dualModel.addVar(lb=-1*INF, ub=0)
                v2 = dualModel.addVar(lb=-1*INF, ub=0)
                self.zLessV1[gen, t] = v1
                self.zLessV2[gen, t] = v2
                # linear     v1 = zUp * dual2
                self.addRelaxConstraints(v1, zUp, dual2, BIGM, dualModel)
                # linear    v2 = zDown * dual2
                self.addRelaxConstraints(v2, zDown, dual2, BIGM, dualModel)
                oobbjj.append(-1 * dual2 * gasBase - gasUpper * v1 + gasDown * v2)


        # add uncertainty yu suan ~
        Uncertainty = 4
        dualModel.addConstr(
            sum(self.zUp.flatten()) <= Uncertainty
        )
        dualModel.addConstr(
            sum(self.zDown.flatten()) <= Uncertainty
        )

        assert dualIndexCount == len(dualIndexList)
        dualModel.setObjective( 1*dualModel.getObjective() + 1*sum(oobbjj) )

        if LOG_OUTPUT:
            print('----generate robust Step 3.4 optimize----------------')
        # dualModel.setParam("MIPFocus", 1)
        dualModel.setParam("OutputFlag", 0)
        dualModel.optimize()

        self.zGreat = to_value(self.zUp)
        self.zLess = to_value(self.zDown)


        for aa, bb in zip(self.zGreat.flatten(), self.zLess.flatten()):
            if aa + bb >= 1.01:
                assert 0
        return dualModel.getObjective().getValue()






























    def givenScenario(self, bGreat, bLess):
        # build robust var
        drawModel = gurobi.Model()

        dM = drawModel


        robust_pressure_well       = to_value(self.well_pressure_trans[0])                         # this is known source pressure

        robust_gas_generator_base         = to_value(self.gas_generator)                                                          # this is known generator load
        robust_gas_generator_reserve_up   = to_value(self.gas_generator_reserve_up)              # this is known generator reserve
        robust_gas_generator_reserve_down = to_value(self.gas_generator_reserve_down)            # this is known generator reserve
        robust_gas_generator              = tonp(dM.addVars(self.gas_generator_num, self.T_fre))
        robust_gas_generator_trans        = self.getTranVars(
            robust_gas_generator, time_interval=int(self.T/self.T_fre))

        robust_flow_load                  = to_value(self.gas_load_trans)                                        # this is known load gas flow

        robust_flow_source                = tonp(dM.addVars(self.gas_well_num, self.T))                          # this is unknown source gas flow

        robust_pressure_load              = tonp(dM.addVars(self.gas_load_num, self.T))                        # this is unknown load/gen node pressure
        robust_pressure_generator         = tonp(dM.addVars(self.gas_generator_num, self.T))



        # build network constraints
        MldPsr = np.vstack((
            robust_flow_load           .reshape(-1, 1),         # this is known value
            robust_gas_generator_trans       .reshape(-1, 1),         # this is var
            robust_pressure_well       .reshape(-1, 1)          # this is known value
        ))

        MsrPld = np.vstack((
             robust_flow_source        .reshape(-1, 1),
             robust_pressure_load      .reshape(-1, 1),
             robust_pressure_generator .reshape(-1, 1)       # all these are vars
        ))

        byPassIndex =  robust_flow_source.reshape(-1, 1).shape[0]



        if LOG_OUTPUT:
            print('----generate robust Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count - byPassIndex):
            row = row + byPassIndex
            if row % 100 == 0:
                if LOG_OUTPUT:
                    print('----generate robust Step 1 ' + str(row) + ' ----------------')
            dM.addConstr(
                lhs=MsrPld[row, 0],
                rhs=(self.YY[row, :]).dot(MldPsr)[0],
                sense='=='
            )

        if LOG_OUTPUT:
            print('----generate robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_fre):
                dM.addConstr(
                    lhs = robust_gas_generator[gen, t],
                    rhs = robust_gas_generator_base[gen, int(t/TIME_LONG_TO_FRE)] +
                          bGreat[gen, t] * robust_gas_generator_reserve_up[gen, int(t/TIME_LONG_TO_FRE)] -
                          bLess[gen, t] * robust_gas_generator_reserve_down[gen, int(t/TIME_LONG_TO_FRE)],
                    sense='=='
                )

        if LOG_OUTPUT:
            print('----generate robust Step 3----------------')
        # build pressure limit
        # node_collection = robust_pressure_load.flatten().tolist() + \
        #                   robust_pressure_generator.flatten().tolist()
        # for load in node_collection:
        #     dM.addConstr(
        #         load <= self.gas_node_pressure_max
        #     )
        #     dM.addConstr(
        #         load >= self.gas_node_pressure_min
        #     )


        if LOG_OUTPUT:
            print('----generate robust Step 3.4 optimize----------------')
        # dualModel.setParam("MIPFocus", 1)
        dM.setParam("OutputFlag", 1)
        dM.optimize()

        rpg = to_value(robust_pressure_generator)
        # dddd = rpg[
        #     rpg[0, 200],
        #     rpg[0, 400],
        #     rpg[0, 600],
        #     rpg[0, 800],
        #     rpg[0, 1000],
        #     rpg[0, 1200-1]
        # ]
        dddd = rpg
        gDDD.append(rpg)

        def drawww():
            import matplotlib.pyplot as plt
            import matplotlib.pyplot as plt
            font1 = {'family': 'Times New Roman',
                     'weight': 'normal',
                     'size': 18,
                     }
            font = {'size': 16}
            # fig = plt.figure(figsize=(8,6))
            # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签

            # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
            fig = plt.figure(figsize=(6, 5))
            # plt.figure()
            ax = fig.add_subplot()
            plt.tick_params(labelsize=12)
            # plt.style.use('ieee')
            plt.plot(dddd.flatten(), '--b', linewidth = 3)
            plt.xlabel('tau', font1)
            plt.ylabel('node pressure (10^(-3) Mpa)', font1)
            plt.title('gas generator node pressure ', font1)
            x=np.linspace(0,1200,5)
            y=80 * np.ones(5)
            plt.plot(x,y, color='red', linewidth=2)
            # plt.text(0,0,'gas pressure lower limit', fontsize=10)
            ax.text(0,81,'gas pressure lower limit', font1)
            plt.grid()
            plt.show()

        drawww()
        return








    #
    def addAppendVars(self):
        self.append_pressure_well         = self.well_pressure_trans
        self.append_gas_generator         = self.model.addVars(self.gas_generator_num, self.T_fre)
        self.append_gas_generator_trans   = self.getTranVars(self.append_gas_generator, time_interval=int(self.T/self.T_fre))

        self.append_flow_source           = tonp(self.model.addVars(self.gas_well_num, self.T))
        self.append_pressure_load         = tonp(self.model.addVars(self.gas_load_num, self.T))
        self.append_pressure_generator    = tonp(self.model.addVars(self.gas_generator_num, self.T))
    def appendConstraintNetwork(self):
        # build network constraints
        self.addAppendVars()
        MldPsr = np.vstack((
            self.gas_load_trans             .reshape(-1, 1),         # this is known value
            self.append_gas_generator_trans .reshape(-1, 1),         # this is var
            self.append_pressure_well       .reshape(-1, 1)          # this is known value
        ))

        MsrPld = np.vstack((
            self.append_flow_source        .reshape(-1, 1),
            self.append_pressure_load      .reshape(-1, 1),
            self.append_pressure_generator .reshape(-1, 1)       # all these are vars
        ))
        byPassIndex = self.append_flow_source.reshape(-1, 1).shape[0]

        if LOG_OUTPUT:
            print('----append constr Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count - byPassIndex):
            row = row + byPassIndex
            if row % 100 == 0:
                if LOG_OUTPUT:
                    print('------- append constraint' + str(row) + '----------')
            self.model.addConstr(
                MsrPld[row, 0] ==
                (self.YY[row, :]).dot(MldPsr)[0]
            )

        if LOG_OUTPUT:
            print('----append constr robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_fre):
                self.model.addConstr(
                    self.append_gas_generator[gen, t] ==
                    (self.gas_generator[gen, int(t/TIME_LONG_TO_FRE)] +
                    self.zGreat[gen, t] * self.gas_generator_reserve_up[gen, int(t/TIME_LONG_TO_FRE)] -
                    self.zLess[gen, t] * self.gas_generator_reserve_down[gen, int(t/TIME_LONG_TO_FRE)]),
                )

        if LOG_OUTPUT:
            print('----append constr robust Step 3----------------')
        # build pressure limit
        node_collection = self.append_pressure_load.flatten().tolist() + \
                          self.append_pressure_generator.flatten().tolist()
        for load in node_collection:
            self.model.addConstr(
                load <= self.gas_node_pressure_max
            )
            self.model.addConstr(
                load >= self.gas_node_pressure_min
            )












    def drawPressure(self):
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        mpl.rcParams.update(mpl.rcParamsDefault)
        plt.subplot(2, 2, 1)
        plt.plot(to_value(self.gas_pressure_generator).flatten())

        plt.subplot(2, 2, 2)
        plt.step(range(6), to_value(self.gas_generator).flatten())

        plt.subplot(2, 2, 3)
        plt.plot(to_value(self.gas_pressure_load).flatten())

        plt.subplot(2, 2, 4)
        plt.step(range(6), to_value(tonp(self.gas_load)).flatten())

        plt.show()

        return

        n = self.gas_node_num
        row = int(n / 2) if int(n % 2) == 0 else int(n / 2) + 1
        column = 2

        for i in range(n):
            plt.subplot(row, column, i + 1)
            plt.plot(to_value(self.node_pressure_trans)[i])
        plt.title('Iter' + str(self.picture_index))
        self.picture_index  = self.picture_index + 1
        plt.show()
    def testMatrix(self):
        import matplotlib.pyplot as plt

        self.Matrix_ab[np.abs(self.Matrix_ab) <= 1e-5 ] = 0
        self.Matrix_cd[np.abs(self.Matrix_cd) <= 1e-5 ] = 0

        self.Matrix_ab = np.around(self.Matrix_ab, 5)
        self.Matrix_cd = np.around(self.Matrix_cd, 5)


        ab = self.Matrix_ab[:, self.node_counts - 1, :]
        cd = self.Matrix_cd[:, self.node_counts - 1, :]
        shape = ab.shape

        time_length = int(shape[0])

        flow2 = np.hstack((np.ones(int(time_length/3)) * 1.5,
                          3 * np.ones(time_length - int(time_length / 3))))
        pressure1 = np.ones(time_length) * 100

        f2p1 = np.hstack((flow2, pressure1))

        f1 = np.dot(ab, f2p1)
        p2 = np.dot(cd, f2p1)

        plt.subplot(1, 2, 1)
        plt.plot(f1)
        plt.subplot(1, 2, 2)
        plt.plot(p2)
        plt.show()


class GenDual:
    def __init__(self, originalModel, addOrigModel=False):
        self.dualModel = gurobi.Model()
        self.origModel = originalModel
        self.exprs = []
        self.senses = []
        self.dualVars = []
        self.origVars = []
        self.consCount = 0
        self.i = []
        self.j = []
        self.v = []
        self.valueLeft = []
        self.coeffMatrix = None
        self.robustIndex = []
        self.addOrigModel = addOrigModel
        self.isSOC = {}
        # self.robustWindOutput = tonp(self.dualModel.addVars(self.origModel.gas_generator_num, self.origModel.T))

    def preset_var_index(self, var_tuple: tuple):
        for var_list in var_tuple:
            self.origVars.extend(var_list)

    def addConstr(self, exprLeft, exprRight, sense, name='XXXXXXXXX', appendIndex=False,
                  coeffInput=None, varInput=None, constantInput=None, varIndex=None, showTime=False):
        if sense == '<=':
            expr = exprLeft - exprRight
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)

        elif sense == '>=':
            expr = exprRight - exprLeft
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)

        elif sense == '==':
            expr = exprLeft - exprRight
            self.addConstr(expr, 0, '>=', 'Great'+name, appendIndex, coeffInput, varInput, constantInput, showTime)
            self.addConstr(expr, 0, '<=', 'Less'+name, appendIndex, coeffInput, varInput, constantInput, showTime)
            return
        else:
            assert 0


        expr = -1 * expr  # guarantee expr >= 0
        t1 = time.time()

        # one dual var per constraint
        dual_var = self.dualModel.addVar(lb=-1 * gurobi.GRB.INFINITY, ub=0, name='Dual'+name)
        self.dualVars.append(dual_var)

        if coeffInput is None:
            # decompose expr as (vars, coeffs, constant) WITHOUT zero-coeff
            cons_vars, coeff, constant = self.getVarsDiff(expr)
            # 断言 vars 里面 没有 相同的 变量，即 合并同类项 的 结果:
            self.assertDiffVars(cons_vars)
            # we maintain a original var order list
            # 扩展 变量 列表，同时 返回 cons_vars 对应的 序列号
            colIndexRet = self.extendVars(self.origVars, cons_vars)
        else:                          # else we input all thing
            assert (len(coeffInput) == len(varIndex))
            cons_vars = varInput
            coeff = coeffInput
            constant = constantInput
            colIndexRet = varIndex

        if showTime:
            t2 = time.time()
            print('stage1 ====>' + str(t2 - t1))

        # [CONSTRUCT MATRIX] get row index, eg. the nth constraint
        rowIndex = self.consCount
        self.consCount = self.consCount + 1

        # [CONSTRUCT MATRIX] get the column index, eg, the var index
        # colIndex = self.getVarIndex(cons_vars)

        # [CONSTRUCT MATRIX] construct the sparse one row for the matrix
        self.i.extend([rowIndex] * len(colIndexRet))
        self.j.extend(colIndexRet)
        self.v.extend(coeff)

        # TODO : SO this is Ax + b >= 0 ?? i forget!!!
        # [CONSTRUCT MATRIX] construct the right hand side
        self.valueLeft.append(constant * 1)

        if appendIndex:
            self.robustIndex.append(rowIndex)

    # pf^2 + qf^2 <= VI
    def addConstrQ(self, plusLeft1, plusLeft2, exprRight, sense):
        if sense == '<=':
            expr = plusLeft1 + plusLeft2 - exprRight
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)
        else:
            assert 0

        var1 = plusLeft1.getVar1(0)
        var2 = plusLeft2.getVar1(0)
        var3 = exprRight.getVar1(0)
        var4 = exprRight.getVar2(0)

        dual1 = self.dualModel.addVar()
        dual2 = self.dualModel.addVar()
        dual3 = self.dualModel.addVar()
        dual4 = self.dualModel.addVar()

        self.isSOC[var1] = True
        self.isSOC[var2] = True
        self.isSOC[var3] = True
        self.isSOC[var4] = True

        self.dualModel.addConstr(
            dual1 * dual1 + dual2 * dual2 <= dual3 * dual4
        )

        # [CONSTRUCT MATRIX] get row index, eg. the nth constraint
        rowIndex = self.consCount
        self.consCount = self.consCount + 1

    def assertDiffVars(self, varList):
        for var1 in varList:
            count = 0
            for var2 in varList:
                if var1.sameAs(var2):
                    count = count + 1
            assert count <= 1

    def getDual(self):
        if LOG_OUTPUT:
            print('----generate robust Step 3.2.1 generate full matrix---------------')
        from scipy.sparse import csr_matrix
        sp = csr_matrix((self.v, (self.i, self.j)), shape=(self.consCount, len(self.origVars)))

        self.coeffMatrix = np.array(sp.toarray())
        self.valueArray = np.array(self.valueLeft)

        self.origVar = np.array(self.origVars)
        self.dualVar = np.array(self.dualVars)

        self.coeffMatrixTrans = self.coeffMatrix.transpose()

        self.objCoeffMatrix = self.objCoeffMatrix

        if LOG_OUTPUT:
            print('----generate robust Step 3.2.2 add constraints----------------')
        for i in range(len(self.origVars)):
            # write a sparse form
            coeff = self.coeffMatrixTrans[i]
            non_zero_index = np.where(coeff != 0)

            coeff_left = coeff[non_zero_index]
            var_right = self.dualVar[non_zero_index]
            # dual_expr_left = sum(self.coeffMatrixTrans[i] * self.dualVar)
            dual_expr_left = sum(coeff_left * var_right)
            dual_expr_right = self.objCoeffMatrix[i]
            if coeff_left.size:
                self.dualModel.addConstr(dual_expr_left == -1*dual_expr_right)
            # if self.senses[i] == '==':
            #     self.dualModel.addConstr(dual_expr_left <= dual_expr_right)
            # if self.senses[i] == '<=':
            #     self.dualModel.addConstr(dual_expr_left <= dual_expr_right)

        self.dualModel.setObjective(-1 * sum(self.valueArray * self.dualVar))

        if LOG_OUTPUT:
            print('----generate robust Step 3.2.3 update model----------------')
        self.dualModel.update()
        return self.dualModel

    def getVars(self, expression):
        assert (isinstance(expression, gurobi.LinExpr))
        self.doNothing = 1

        expr = expression
        varSize = expr.size()
        allVars = []
        allCoeff = []

        for index in range(varSize):
            allVars.append(expr.getVar(index))
            allCoeff.append(expr.getCoeff(index))

        constant = expr.getConstant()
        return allVars, allCoeff, constant

    def getVarsDiff(self, expression):
        assert (isinstance(expression, gurobi.LinExpr))

        expr = expression
        varSize = expr.size()
        allVars = []
        allCoeff = []

        for index in range(varSize):
            allVars.append(expr.getVar(index))
            allCoeff.append(expr.getCoeff(index))

        constant = expr.getConstant()

        var_value = {}
        self.origModel.update()
        for index, var in enumerate(allVars):
            if var not in var_value:
                var_value[var] = allCoeff[index]
            else:
                var_value[var] += allCoeff[index]

        allVars = list(var_value.keys())
        allCoeff = list(var_value.values())
        return allVars, allCoeff, constant


    def extendVars(self, listTo, all_vars):
        insertIndex = []
        for var in all_vars:
            index = self.getIndex(listTo, var)
            if index == -1:
                listTo.append(var)
                insertIndex.append(len(listTo) - 1)
            else:
                insertIndex.append(index)
        return insertIndex

    def getVarIndex(self, cons_vars):
        # self.origModel.update()
        indices = []
        for var in cons_vars:
            pos = self.getIndex(self.origVars, var)
            indices.append(pos)
        return indices

    def getIndex(self, varList, var):
        self.doNothing = 1
        for ind, varAt in enumerate(varList):
            if varAt.sameAs(var):
                return ind
        return -1

    def addObjectiveMin(self, expr):
        from scipy.sparse import csr_matrix
        obj_var_len = expr.size()
        obj_i = []
        obj_j = []
        obj_v = []

        for index in range(obj_var_len):
            var = expr.getVar(index)
            ind = self.getIndex(self.origVars, var)
            coe = expr.getCoeff(index)
            if ind == -1:
                assert 0

            obj_i.append(0)
            obj_j.append(ind)
            obj_v.append(coe)

        sp = csr_matrix((obj_v, (obj_i, obj_j)), shape=(1, len(self.origVars)))

        self.objCoeffMatrix = np.array(sp.toarray())[0]

        if self.addOrigModel:
            self.origModel.setObjective(1 * expr)

class GenDual2:
    def __init__(self, originalModel, addOrigModel=False):
        self.dualModel = gurobi.Model()
        self.origModel = originalModel
        self.exprs = []
        self.senses = []
        self.dualVars = []
        self.origVars = []
        self.consCount = 0
        self.i = []
        self.j = []
        self.v = []
        self.valueLeft = []
        self.coeffMatrix = None
        self.robustIndex = []
        self.addOrigModel = addOrigModel
        self.isSOC = {}
        self.socMap = {}

    def preset_var_index(self, var_tuple: tuple):
        for var_list in var_tuple:
            self.origVars.extend(var_list)

    def addConstr(self, exprLeft, exprRight, sense, name='XXXXXXXXX', appendIndex=False,
                   ):
        if sense == '<=':
            expr = exprLeft - exprRight
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)

        elif sense == '>=':
            expr = exprRight - exprLeft
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)

        elif sense == '==':
            expr = exprLeft - exprRight
            self.addConstr(expr, 0, '>=', 'Great' + name, appendIndex,  )
            self.addConstr(expr, 0, '<=', 'Less' + name, appendIndex,  )
            return
        else:
            assert 0

        expr = -1 * expr  # guarantee expr >= 0

        # one dual var per constraint
        dual_var = self.dualModel.addVar(lb=-1 * gurobi.GRB.INFINITY, ub=0, name='Dual' + name)
        self.dualVars.append(dual_var)

        # decompose expr as (vars, coeffs, constant) WITHOUT zero-coeff
        cons_vars, coeff, constant = self.getVarsDiff(expr)
        # 断言 vars 里面 没有 相同的 变量，即 合并同类项 的 结果:
        self.assertDiffVars(cons_vars)
        # we maintain a original var order list
        # 扩展 变量 列表，同时 返回 cons_vars 对应的 序列号
        colIndexRet = self.extendVars(self.origVars, cons_vars)


        # [CONSTRUCT MATRIX] get row index, eg. the nth constraint
        rowIndex = self.consCount
        self.consCount = self.consCount + 1

        # [CONSTRUCT MATRIX] get the column index, eg, the var index
        # colIndex = self.getVarIndex(cons_vars)

        # [CONSTRUCT MATRIX] construct the sparse one row for the matrix
        self.i.extend([rowIndex] * len(colIndexRet))
        self.j.extend(colIndexRet)
        self.v.extend(coeff)

        # TODO : SO this is Ax + b >= 0 ?? i forget!!!
        # [CONSTRUCT MATRIX] construct the right hand side
        self.valueLeft.append(constant * 1)

        if appendIndex:
            self.robustIndex.append(rowIndex)

    # pf^2 + qf^2 <= VI
    def addConstrQ(self, plusLeft1, plusLeft2, plusLeft3, exprRight, sense):
        if sense == '<=':
            expr = plusLeft1 + plusLeft2 + plusLeft3 - exprRight
            self.exprs.append(-1 * expr)
            self.senses.append('>=')
            if self.addOrigModel:
                self.origModel.addConstr(expr <= 0)
        else:
            assert 0

        var1 = plusLeft1.getVar1(0)   # a
        var2 = plusLeft2.getVar1(0)   # b
        var3 = plusLeft3.getVar1(0)   # c
        var4 = exprRight.getVar2(0)   # d   a^2 + b^2 + c^2 <= d^2

        dual1 = self.dualModel.addVar(lb=-1*INF, ub=INF)
        dual2 = self.dualModel.addVar(lb=-1*INF, ub=INF)
        dual3 = self.dualModel.addVar(lb=-1*INF, ub=INF)
        dual4 = self.dualModel.addVar()

        self.isSOC[var1] = True
        self.isSOC[var2] = True
        self.isSOC[var3] = True
        self.isSOC[var4] = True

        self.socMap[var1] = dual1
        self.socMap[var2] = dual2
        self.socMap[var3] = dual3
        self.socMap[var4] = dual4

        self.dualModel.addConstr(
            dual1 * dual1 + dual2 * dual2 + dual3 * dual3 <= dual4 * dual4
        )

        self.extendVars(self.origVars, [var1, var2, var3, var4])
        # instead extend dualVars, we just store dual var to map and add linear constrain later


    def getDual(self):
        if LOG_OUTPUT:
            print('----generate robust Step 3.2.1 generate full matrix---------------')
        from scipy.sparse import csr_matrix
        sp = csr_matrix((self.v, (self.i, self.j)), shape=(self.consCount, len(self.origVars)))

        self.coeffMatrix = np.array(sp.toarray())
        self.valueArray = np.array(self.valueLeft)

        self.origVar = np.array(self.origVars)
        self.dualVar = np.array(self.dualVars)

        self.coeffMatrixTrans = self.coeffMatrix.transpose()

        self.objCoeffMatrix = self.objCoeffMatrix

        if LOG_OUTPUT:
            print('----generate robust Step 3.2.2 add constraints----------------')

        for i in range(len(self.origVars)):
            if not self.origVars[i] in self.isSOC:
                # write a sparse form
                coeff = self.coeffMatrixTrans[i]
                non_zero_index = np.where(coeff != 0)

                coeff_left = coeff[non_zero_index]
                var_right = self.dualVar[non_zero_index]
                # dual_expr_left = sum(self.coeffMatrixTrans[i] * self.dualVar)
                dual_expr_left = sum(coeff_left * var_right)
                dual_expr_right = self.objCoeffMatrix[i]
                if coeff_left.size:
                    self.dualModel.addConstr(dual_expr_left == -1 * dual_expr_right)
            else:
                coeff = self.coeffMatrixTrans[i]
                non_zero_index = np.where(coeff != 0)
                coeff_left = coeff[non_zero_index]
                var_right = self.dualVar[non_zero_index]
                dual_expr_left = sum(coeff_left * var_right)
                dual_expr_right = self.objCoeffMatrix[i]
                dual_expr = dual_expr_left + dual_expr_right
                self.dualModel.addConstr(dual_expr == 1 * self.socMap[self.origVars[i]])
        self.dualModel.setObjective(-1 * sum(self.valueArray * self.dualVar))

        if LOG_OUTPUT:
            print('----generate robust Step 3.2.3 update model----------------')
        self.dualModel.update()
        return self.dualModel




    def assertDiffVars(self, varList):
        for var1 in varList:
            count = 0
            for var2 in varList:
                if var1.sameAs(var2):
                    count = count + 1
            assert count <= 1


    def getVarsDiff(self, expression):
        assert (isinstance(expression, gurobi.LinExpr))

        expr = expression
        varSize = expr.size()
        allVars = []
        allCoeff = []

        for index in range(varSize):
            allVars.append(expr.getVar(index))
            allCoeff.append(expr.getCoeff(index))

        constant = expr.getConstant()

        var_value = {}
        self.origModel.update()
        for index, var in enumerate(allVars):
            if var not in var_value:
                var_value[var] = allCoeff[index]
            else:
                var_value[var] += allCoeff[index]

        allVars = list(var_value.keys())
        allCoeff = list(var_value.values())
        return allVars, allCoeff, constant

    def extendVars(self, listTo, all_vars):
        insertIndex = []
        for var in all_vars:
            index = self.getIndex(listTo, var)
            if index == -1:
                listTo.append(var)
                insertIndex.append(len(listTo) - 1)
            else:
                insertIndex.append(index)
        return insertIndex

    def getIndex(self, varList, var):
        for ind, varAt in enumerate(varList):
            if varAt.sameAs(var):
                return ind
        return -1

    def addObjectiveMin(self, expr):
        from scipy.sparse import csr_matrix
        obj_var_len = expr.size()
        obj_i = []
        obj_j = []
        obj_v = []

        for index in range(obj_var_len):
            var = expr.getVar(index)
            ind = self.getIndex(self.origVars, var)
            coe = expr.getCoeff(index)
            if ind == -1:
                assert 0

            obj_i.append(0)
            obj_j.append(ind)
            obj_v.append(coe)

        sp = csr_matrix((obj_v, (obj_i, obj_j)), shape=(1, len(self.origVars)))

        self.objCoeffMatrix = np.array(sp.toarray())[0]

        if self.addOrigModel:
            self.origModel.setObjective(1 * expr)

def testGetDual():
    model = gurobi.Model()
    x1 = model.addVar(name='x1')
    x2 = model.addVar(name='x2')
    x3 = model.addVar(name='x3')
    x4 = model.addVar(name='x4')
    x5 = model.addVar(name='x5')

    gD = GenDual(model, addOrigModel=True)
    gD.addConstr(42.4 * x1 + 3 * x2, 90, '<=', name='c1')
    gD.addConstr(2.3 * x1 + x2, 80, '<=', name='c2')
    gD.addConstr(2.2 * x3 + 33*x2, 80, '==', name='c3')
    gD.addConstr(x1 + 4*x2, 25, '<=', name='c4')
    gD.addConstr(x1, 0, '>=', name='c5')
    gD.addConstr(x2, 0, '>=', name='c6')
    gD.addConstr(x3, 0, '>=', name='c7')
    #
    gD.addObjectiveMin(-3.5 * x1 - 4 * x2)
################
    # gD.addConstr(x1 + x2 + 2 * x3 + x4 + 3 * x5, 4, '>=')
    # gD.addConstr(2 * x1 - x2 + 3 * x3 + x4 + x5, 3, '>=')
    # gD.addConstr(x1 + x4, 4, '==')
    # gD.addConstr(x1, 0, '>=')
    # gD.addConstr(x2, 0, '>=')
    # gD.addConstr(x3, 0, '>=')
    # gD.addConstr(x4, 0, '>=')
    # gD.addConstr(x5, 0, '>=')
    # gD.addObjectiveMin((2 * x1 + 3 * x2 + 5 * x3 + 2 * x4 + 3 * x5) * 1)
#################
    # gD.addConstr(x1 + 2 * x2, 2, '<=')
    # gD.addConstr(x1 - x2, 3, '<=')
    # gD.addConstr(2 * x1 + 3 * x2, 5, '<=')
    # gD.addConstr(x1 + x2, 2, '<=')
    # gD.addConstr(3 * x1 + x2, 3, '<=')
    # gD.addConstr(x1, 0, '>=')
    # gD.addConstr(x2, 0, '>=')
    # gD.addObjectiveMin(-4 * x1 - 3 * x2)

    dual = gD.getDual()
    dual.update()

    model.write('originalModel.lp')
    dual.write('dualModel.lp')
    dual.optimize()
    model.optimize()

    print(dual.getObjective().getValue() )
    print(model.getObjective().getValue() * -1)


def doMain():
    pg.buildBasePrimaryProblemNetwork()

    # pg.model.optimize()
    # pg.feasibleProblem()
    # pg.appendConstraint()


    for i in range(30):
        print('Stage 1 : preSolve : [ Iteration: ' + str(i) + ']')
        # pg.model.setParam('OutputFlag', 0)
        pg.model.setParam("LogFile", "log" + str(i) + '.txt')
        pg.model.optimize()
        pg.drawPressure()
        print('    ===> Stage 1.1 : objective ' + str(pg.model.getObjective().getValue()))
        print('    ===> Stage 2 : feasibleTest : Reserve Up [' + str(to_value(pg.gas_generator_reserve_up)) + '] Reserve Down [' + str(to_value(pg.gas_generator_reserve_down)) + ']')
        t1 = time.time()
        pg.feasibleProblemNetworkPreset()
        feasibleTest = pg.feasibleProblemNetworkBuild(pg.robustModel)
        # pg.testOriginalFeasibelModel()
        # feasibleTest = pg.feasibleProblemNetworkTest()
        gFeasibleTest.append(feasibleTest)
        gProfile.append(pg.model.getObjective().getValue())
        t2 = time.time()
        if LOG_OUTPUT:
            print('Feasibel Time : ' + str(t2 - t1))
        if abs(feasibleTest) <= 3e2:
            print('    ===> Stage 3 : Terminate : feasible test_pg result : [' + str(feasibleTest) + ']')
            print('    ===> feasibleTest : Reserve Up [' + str(
                to_value(pg.gas_generator_reserve_up)) + '] Reserve Down [' + str(
                to_value(pg.gas_generator_reserve_down)) + ']')

            return
        else:
            print('    ===> Stage 3 : Append constraints : feasible test_pg result : [' + str(feasibleTest) + ']')
            print('    ===> Stage 3 : Scenario found is : ')
            print('           ' + str(pg.zLess.flatten().tolist()))
            print('           ' + str(pg.zGreat.flatten().tolist()))
            scenarioLess.append(pg.zLess.tolist())
            scenarioGreat.append(pg.zGreat.tolist())
            pg.givenScenario(pg.zGreat, pg.zLess)
            pg.appendConstraintNetwork()
        print(' ')



def doTest():
    # 1, test_pg network form generation correction
    from test_pg.network6node import getConfig as testConfig
    pg_test = PowerGas(*testConfig())





def drawScenario():
    import matplotlib.pyplot as plt
    plt.figure(figsize=(8,8))

    sLen = len(scenarioLess)
    row = int(sLen / 2) if sLen % 2 == 0 else int(sLen / 2 + 1)

    for i in range(sLen):
        plt.subplot(row, 2, i+1)
        plt.plot(scenarioGreat[i])
        plt.plot(scenarioLess[i])

    plt.show()

def drawSce():
    import matplotlib.pyplot as plt
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 14,
             }
    plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签

    plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    plt.figure(figsize=(8,4))
    # plt.style.use( 'ieee')

    sLen = 2# len(scenarioLess)
    row = int(sLen / 2) if sLen % 2 == 0 else int(sLen / 2 + 1)

    for i in range(sLen):
        plt.subplot(row, 2, i+1)
        plt.step([0,1,2,3,4,5,6,7,8,9,10,11], scenarioGreat[i][0], color='black', linewidth=3, label='upper')
        plt.step([0,1,2,3,4,5,6,7,8,9,10,11], scenarioLess[i][0], '--', color='red', linewidth= 3,label='down')
        plt.xlabel('time slot', font1)
        plt.ylabel('reserve demand', font1)
        plt.title('PD from power grid', font1)
        plt.legend()

    plt.show()

def printIJV(row, dG):
    index = np.where(np.array(dG.i) == row)
    jj = np.array(dG.j)[index]
    vv = np.array(dG.v)[index]
    print(jj)
    print(vv)

def drawReserveAndBaseValue():
    import matplotlib.pyplot as plt
    import matplotlib as mpl
    plt.tight_layout()

    mpl.rcParams.update(mpl.rcParamsDefault)
    plt.figure(figsize=(5, 4.4))
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 18,
             }
    # plt.style.use( 'ieee')

    rUp = to_value(pg.gas_generator_reserve_up)
    rDown = to_value(pg.gas_generator_reserve_down)
    base = to_value(pg.gas_generator)

    cap = pg.gas_gen_capacity * pg.T_long
    if rUp.shape[0] == 1:
        plt.plot(cap, label='capacity', linewidth = 2, color=
                 'darkviolet')
        plt.plot(rUp[0], ':', label='upper reserve', linewidth = 2 , color='maroon')
        plt.plot([0]*6, '-.', label='lower reserve' , linewidth = 2)
        plt.plot(base[0], '--*', label='operation point' , linewidth = 3, color='black', markersize=10)

        plt.xlabel('time', font1 )
        plt.ylabel('power output (MW)',font1 )
        # plt.title('Gas generator reserve and operation point')
        plt.legend(loc='upper right', prop={'size':14} )
        plt.title('gas generator operation point and reserve ', font1)
        plt.grid()
        plt.show()
    else:
         assert 0


def drawAlgorithm():
    import matplotlib.pyplot as plt
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 18,
             }
    # plt.style.use('ieee')
    plt.figure(figsize=(7, 5))
    aaa = plt.plot([-1 * data for data in gFeasibleTest], '--*', linewidth=3, markersize=15, label='feasible result', color='red')
    # plt.plot(gProfile, label='IES profile')

    plt.xlabel('iteration index', font1)
    plt.ylabel('feasible result', font1)
    plt.title('robust feasibility / objective dynamic', font1)
    # plt.legend()
    plt.grid()
    plt.twinx()
    gObj = [33.2, 30.3, 29.1, 29]
    bbb = plt.plot([-1*f for f in gProfile], '--.', linewidth=3, markersize=15, label='objective', color='blue')
    plt.ylim((5,25))
    plt.ylabel('Profits of IEGS (*10^4 $)', font1)
    labss = [l.get_label() for l in aaa + bbb]
    plt.legend(aaa + bbb, labss, loc=0, prop={'size':16} )
    # plt.legend()
    plt.show()


def drawPressureDynamic():
    plt.figure(figsize=(7,6))
    plt.tight_layout()
    plt.subplots_adjust(wspace=.3, hspace=.5)
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 14,
             }
    def getXXX(a) :
        a1 = a[0]
        return [a1] + a

    plt.subplot(2, 2, 1)
    plt.plot(to_value(pg.gas_pressure_generator).flatten(), color='black', linewidth = 2)
    plt.xlabel('tau (0.75s)', font1)
    plt.ylabel('node pressure (10^(-3) Mpa)', font1)
    plt.title('gas generator pressure', font1)
    plt.grid()
    plt.subplot(2, 2, 2)
    plt.step(range(7), getXXX(to_value(pg.gas_generator).flatten().tolist()), color='black', linewidth = 2)
    plt.xlabel('t (15min)', font1)
    plt.ylabel('gas flow (Kg/s)', font1)
    plt.title('gas generator flow', font1)
    plt.grid()

    plt.subplot(2, 2, 3)
    plt.plot(to_value(pg.gas_pressure_load).flatten(), color='black', linewidth = 2)
    plt.xlabel('tau (0.75s)', font1)
    plt.ylabel('node pressure (10^(-3) Mpa)', font1)
    plt.title('gas load pressure', font1)
    plt.grid()

    plt.subplot(2, 2, 4)

    plt.step(range(7), getXXX(to_value(tonp(pg.gas_load)).flatten().tolist()), color='black', linewidth = 2)
    plt.xlabel('t (15min)', font1)
    plt.ylabel('gas flow (Kg/s)', font1)
    plt.title('gas load flow', font1)
    plt.grid()

    plt.show()

def drawSecAndDynamic():
    import matplotlib.pyplot as plt
    font1 = {'family': 'Times New Roman',
             'weight': 'normal',
             'size': 16,
             }

    # plt.rcParams['font.sans-serif'] = ['SimHei']  # 用来正常显示中文标签
    #
    # plt.rcParams['axes.unicode_minus'] = False  # 用来正常显示负号
    plt.figure(figsize=(8,7))
    # plt.style.use( 'ieee')
    plt.tight_layout()
    plt.subplots_adjust(wspace=.3, hspace=.3)

    plt.subplot(2, 2, 1)
    plt.step((np.array([0,1,2,3,4,5,6,7,8,9,10,11])/2).tolist(), np.roll(np.array(scenarioGreat[0][0]), 1), color='black', linewidth=3, label='upper')
    plt.step((np.array([0,1,2,3,4,5,6,7,8,9,10,11])/2).tolist(), np.roll(np.array(scenarioLess[0][0]), 1), '--', color='red', linewidth= 3,label='down')
    plt.xlabel('time slot', font1)
    plt.ylabel('reserve demand', font1)
    plt.title('a. PD from power grid', font1)
    plt.legend()
    plt.grid()

    plt.subplot(2, 2, 2)
    plt.step((np.array([0,1,2,3,4,5,6,7,8,9,10,11])/2).tolist(), np.roll(np.array(scenarioGreat[1][0]), 1), color='black', linewidth=3, label='upper')
    plt.step((np.array([0,1,2,3,4,5,6,7,8,9,10,11])/2).tolist(), np.roll(np.array(scenarioLess[1][0]), 1), '--', color='red', linewidth= 3,label='down')
    plt.xlabel('time slot', font1)
    plt.ylabel('reserve demand', font1)
    plt.title('b. PD from power grid', font1)
    plt.legend()
    plt.grid()




    ax = plt.subplot(2, 2, 3)
    dddd = gDDD[0]
    plt.tick_params(labelsize=12)
    # plt.style.use('ieee')
    plt.plot(dddd.flatten(), '--b', linewidth=3)
    plt.xlabel('tau', font1)
    plt.ylabel('node pressure (10^(-3) Mpa)', font1)
    plt.title('c. gas generator node pressure ', font1)
    x = np.linspace(0, 1200, 5)
    y = 80 * np.ones(5)
    plt.plot(x, y, color='red', linewidth=2)
    # plt.text(0,0,'gas pressure lower limit', fontsize=10)
    ax.text(0, 81, 'gas pressure lower limit', font1)
    plt.grid()





    ax = plt.subplot(2, 2, 4)
    dddd = gDDD[1]
    plt.tick_params(labelsize=12)
    # plt.style.use('ieee')
    plt.plot(dddd.flatten(), '--b', linewidth=3)
    plt.xlabel('tau', font1)
    plt.ylabel('node pressure (10^(-3) Mpa)', font1)
    plt.title('d. gas generator node pressure ', font1)
    x = np.linspace(0, 1200, 5)
    y = 80 * np.ones(5)
    plt.plot(x, y, color='red', linewidth=2)
    # plt.text(0,0,'gas pressure lower limit', fontsize=10)
    ax.text(0, 81, 'gas pressure lower limit', font1)
    plt.grid()

    plt.show()

if __name__ == '__main__':
    # testGetDual()
    # assert 0
    gFeasibleTest = []
    gProfile = []
    gDDD = []
    pg = PowerGas(*getConfig())

    # pg.testMatrix()
    # assert 0


    scenarioLess = []
    scenarioGreat = []

    doMain()
    # doTest()