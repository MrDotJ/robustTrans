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

LOG_OUTPUT = False
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





    def build_original_model_var(self):
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

        self.power_gen_trans = self.getTranVars(self.power_gen)
        self.react_gen_trans = self.getTranVars(self.react_gen)
        self.power_load_trans = self.getTranVars(self.power_load)
        self.react_load_trans = self.getTranVars(self.react_load)

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
        # 线路流量
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
        # self.start_pipeline_pressure_trans = tonp(self.model.addVars(
        #     self.gas_line_num, self.T,
        #     lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
        #     name='pipeline_gas_pressure_of_start_port'
        # ))
        # self.end_pipeline_pressure_trans = tonp(self.model.addVars(
        #     self.gas_line_num, self.T,
        #     lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
        #     name='pipeline_gas_pressure_of_end_port'
        # ))

        ########
        # 节点气压
        self.node_pressure_trans = tonp(self.model.addVars(
            self.gas_node_num, self.T,
            lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
            name='gas_node_pressure_trans'
        ))
        ########

        # output of gas generator is T-long, but the reserve is real-time
        # 燃气轮机的 出力 以及 备用， 小时级
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
        # self.well_pressure = self.model.addVars(
        #     self.gas_well_num, self.T_long,
        #     name='gas_well_set_pressure'
        # )
        self.gas_generator_trans = self.getTranVars(self.gas_generator)
        self.gas_load_trans = self.getTranVars(self.gas_load)
        # self.well_pressure_trans = getTranVars(self.well_pressure)
    def buildBasePrimaryProblem(self):
        # ----------- 1. add Vars ---------------------------
        self.build_original_model_var()
        # ----------- 2.1 node power balance -----------------
        self.model.update()
        #
        # print('start')
        # for node in range(self.bus_num):
        #     for time in range(self.T_long):
        #         self.model.addConstr(
        #             sum(self.power_gen[np.where(np.array(self.gen_conn_power_index) == node), time].flatten()) +
        #             sum(self.line_power_flow[np.where(np.array(self.line_end_point) == node), time].flatten()) -
        #             sum((self.line_resistance[np.where(np.array(self.line_end_point) == node)] *
        #                  self.line_current_square[np.where(np.array(self.line_end_point) == node), time]).flatten()) +
        #             sum(self.gas_generator[np.where(np.array(self.gas_gen_index_power) == node), time].flatten())
        #             ==
        #             sum(self.power_load[np.where(np.array(self.load_index) == node), time].flatten()) +
        #             sum(self.line_power_flow[np.where(np.array(self.line_start_point) == node), time].flatten()),
        #             name='power_balance')
        #         self.model.addConstr(
        #             sum(self.react_gen[np.where(np.array(self.gen_conn_power_index) == node), time].flatten()) +
        #             sum(self.line_react_flow[np.where(np.array(self.line_end_point) == node), time].flatten()) -
        #             sum((self.line_reactance[np.where(np.array(self.line_end_point) == node)] *
        #                  self.line_current_square[np.where(np.array(self.line_end_point) == node), time]).flatten())
        #             ==
        #             sum(self.react_load[np.where(np.array(self.load_index) == node), time].flatten()) +
        #             sum(self.line_react_flow[np.where(np.array(self.line_start_point) == node), time].flatten()),
        #             name='react_balance')
        #
        # # ----------- 2.2 line voltage drop ------------------
        # for i in range(self.line_num):
        #     start_point, end_point = self.line_start_point[i], self.line_end_point[i]
        #     resistance, reactance = self.line_resistance[i], self.line_reactance[i]
        #     impedance_square = reactance * reactance + resistance * resistance
        #     for time in range(self.T_long):
        #         self.model.addConstr(
        #             self.voltage_square[end_point, time] -
        #             self.voltage_square[start_point, time]
        #             ==
        #             impedance_square * self.line_current_square[i, time] -
        #             2 * (resistance * self.line_power_flow[i, time] +
        #                  reactance * self.line_react_flow[i, time]),
        #             name='voltage_drop')
        #         self.model.addConstr(
        #             self.line_power_flow[i, time] * self.line_power_flow[i, time] +
        #             self.line_react_flow[i, time] * self.line_react_flow[i, time]
        #             <=
        #             self.line_current_square[i, time] * self.voltage_square[start_point, time],
        #             name='flow_relax')

        # ====> Two port network of gas pipeline [gas system transient]

        self.Matrix_ab[np.abs(self.Matrix_ab) <= 1e-5 ] = 0
        self.Matrix_cd[np.abs(self.Matrix_cd) <= 1e-5 ] = 0

        self.Matrix_ab = np.around(self.Matrix_ab, 5)
        self.Matrix_cd = np.around(self.Matrix_cd, 5)

        def myDot2(matrix, syms, row, y_ranges):
            value = matrix[row][self.node_counts - 1][y_ranges]
            sym = syms[y_ranges]
            return sum(value * sym)

        # ====> 2.3 two port pipeline [gas system transient]
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


        #
        #         print(myDot2(self.Matrix_ab, flow2pressure1, t, y_range))
        #         print(myDot2(self.Matrix_cd, flow2pressure1, t, y_range2))
        #
        #
        # assert 0
        # ====> 2.4 node gas balance [gas system transient]
        for t in range(self.T):
            for node in range(self.gas_node_num):
                self.model.addConstr(
                    sum(self.start_pipeline_flow_trans[np.where(self.gas_pipeline_start_node == node),   t].flatten()) +
                    sum(self.gas_load_trans[           np.where(self.gas_load_connection_index == node), t].flatten()) +
                    sum(self.gas_generator_trans[      np.where(self.gas_gen_index_gas == node),         t].flatten())
                    ==
                    sum(self.end_pipeline_flow_trans[  np.where(self.gas_pipeline_end_node == node),     t].flatten()) +
                    sum(self.well_in_trans[            np.where(self.gas_well_connection_index == node), t].flatten()),
                    name='c19')

        # ====> 2.6 gas generator transient [gas system transient]
        # print('======> finish 7.3')
        for t in range(self.T_long):
            for gen in range(self.gas_generator_num):
                self.model.addConstr(
                    self.gas_generator[gen, t] + self.gas_generator_reserve_up[gen, t] <=
                    self.gas_gen_capacity[gen], name='c22'
                )
                self.model.addConstr(
                    self.gas_generator[gen, t] - self.gas_generator_reserve_down[gen, t] >=
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
                objs.append(self.well_in_trans[well, t] * self.well_cost / self.time_counts)

        for load in range(self.gas_load_num):
            for t in range(self.T_long):
                objs.append((self.gas_load[load, t] - (self.gas_load_min[load, t] + self.gas_load_max[load, t]) / 2) *
                            (self.gas_load[load, t] - (self.gas_load_min[load, t] + self.gas_load_max[load, t]) / 2))

        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                objs.append(-1 * self.gas_generator_reserve_up[gen, t] - self.gas_generator_reserve_down[gen, t])
        # --------------- 3 set objective ----------------
        self.model.setParam('OutputFlag', 0)
        self.model.setParam('BarHomogeneous', 1)
        # self.model.setParam('Method', 1)
        self.model.setObjective(sum(objs))






    # [TEST] this function just test_pg the feasibility of original problem which
    # used in feasible problem
    def testOriginalFeasibleModel(self):
        self.testModel = gurobi.Model()
        tM = self.testModel


        self.robust_pressure_well = to_value(self.node_pressure_trans[0])  # this is known source pressure

        self.robust_gas_generator_base = to_value(self.gas_generator) * 1  # this is known generator load
        self.robust_gas_generator_reserve_up = to_value(
            self.gas_generator_reserve_up) * 0  # this is known generator reserve
        self.robust_gas_generator_reserve_down = to_value(
            self.gas_generator_reserve_down) * 0  # this is known generator reserve
        self.robust_gas_generator = tonp(tM.addVars(self.gas_generator_num, self.T_long))
        self.robust_gas_generator_trans = self.getTranVars(self.robust_gas_generator)

        self.robust_flow_load = to_value(self.gas_load_trans)  # this is known load gas flow

        self.robust_flow_source = tonp(tM.addVars(self.gas_well_num, self.T))  # this is unknown source gas flow

        self.robust_pressure_load = tonp(
            tM.addVars(self.gas_load_num, self.T))  # this is unknown load/gen node pressure
        self.robust_pressure_generator = tonp(
            tM.addVars(self.gas_generator_num, self.T))  # this is unknown load/gen node pressure

        self.zUp = tonp(tM.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zDown = tonp(tM.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zUp[0,0] = 1
        self.zUp[0,1] = 0
        self.zDown[0,0] = 0
        self.zDown[0,1] = 1

        MldPsr = np.vstack((
            self.robust_flow_load           .reshape(-1, 1),         # this is known value
            self.robust_gas_generator_trans       .reshape(-1, 1),         # this is var
            self.robust_pressure_well       .reshape(-1, 1)          # this is known value
        ))
        start_index = self.robust_flow_load.reshape(-1, 1).shape[0]
        middle_index = self.robust_flow_load.reshape(-1, 1).shape[0] + self.robust_gas_generator_trans.reshape(-1, 1).shape[0]

        MsrPld = np.vstack((
            self.robust_flow_source        .reshape(-1, 1),
            self.robust_pressure_load      .reshape(-1, 1),
            self.robust_pressure_generator .reshape(-1, 1)       # all these are vars
        ))

        byPassIndex = self.robust_flow_source.reshape(-1, 1).shape[0]

        # because of the O^2 complex of index search , so we pre-set variable index
        MAGIC_LIST = np.arange(self.robust_gas_generator.reshape(-1, 1).shape[0]).tolist()

        map_value_index = 0
        index_m = {}
        self.testModel.update()

        var_preset = np.array(self.robust_gas_generator.reshape(-1, 1).flatten().tolist() +
            self.robust_flow_source.reshape(-1, 1).flatten().tolist() +
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist() +
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist())

        for vars_preset in (
            self.robust_gas_generator.reshape(-1, 1).flatten().tolist(),
            self.robust_flow_source.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist()
        ):
            vars_list = vars_preset
            for var in vars_list:
                index_m[var] = map_value_index
                map_value_index = map_value_index + 1

        print('----generate robust Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count - byPassIndex):
            row = row + byPassIndex
            if row % 100 == 0:
                print('----generate robust Step 1 ' + str(row) + ' ----------------')
            # self.testModel.addConstr(
            #     MsrPld[row, 0] ==
            #     (self.YY[row, :]).dot(MldPsr)[0]
            # )
            # continue
            coeffInput = [ -1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (1 * np.sum(coeff, axis=1)).tolist())

            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex. extend( MAGIC_LIST )

            constantInput = 1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )

            self.testModel.addConstr(
                sum(var_preset[varIndex] * np.array(coeffInput)) + constantInput <= 0
            )

            coeffInput = [ 1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (-1 * np.sum(coeff, axis=1)).tolist())
            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex.extend( MAGIC_LIST )

            constantInput = -1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            self.testModel.addConstr(
                sum(var_preset[varIndex] * np.array(coeffInput)) + constantInput <= 0
            )

        print('----generate robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                self.testModel.addConstr(
                    self.robust_gas_generator[gen, t] == (
                        self.robust_gas_generator_base[gen, t] +
                        self.zUp[gen, t] * self.robust_gas_generator_reserve_up[gen, t] -
                        self.zDown[gen, t] * self.robust_gas_generator_reserve_down[gen, t]
                    ))

        print('----generate robust Step 3----------------')
        # build pressure limit
        sCollection = []
        node_collection = self.robust_pressure_load.flatten().tolist() + \
                          self.robust_pressure_generator.flatten().tolist()
        for load in node_collection:
            sMin = self.testModel.addVar()
            sPlus = self.testModel.addVar()
            sCollection.append(sMin)
            sCollection.append(sPlus)

            self.testModel.addConstr(
                load - sMin <= self.gas_node_pressure_max,
            )
            self.testModel.addConstr(
                load + sPlus >= self.gas_node_pressure_min,
            )

        print('----generate robust Step 3.1 add objective----------------')
        # build primitive objective
        self.testModel.setObjective(gurobi.quicksum(sCollection))

        self.testModel.optimize()
    # [TEST] these two functions test_pg the feasible problem with pre-set scenario
    # that means the zUp and zDown is known value
    # and we don't need the additional objective
    def feasibleProblemNetworkTest(self):
        print('test_pg ==== begin')
        # preset set variable
        self.testModel = gurobi.Model()
        self.dualGen = GenDual(self.testModel, addOrigModel=True)
        tM = self.testModel

        self.robust_pressure_well              = to_value(self.node_pressure_trans[0])                         # this is known source pressure

        self.robust_gas_generator_base         = to_value(self.gas_generator) * 1                                                         # this is known generator load
        self.robust_gas_generator_reserve_up   = to_value(self.gas_generator_reserve_up)  * 0             # this is known generator reserve
        self.robust_gas_generator_reserve_down = to_value(self.gas_generator_reserve_down)   * 0         # this is known generator reserve
        self.robust_gas_generator              = tonp(tM.addVars(self.gas_generator_num, self.T_long))
        self.robust_gas_generator_trans        = self.getTranVars(self.robust_gas_generator)

        self.robust_flow_load                  = to_value(self.gas_load_trans)                                        # this is known load gas flow

        self.robust_flow_source                = tonp(tM.addVars(self.gas_well_num, self.T))                          # this is unknown source gas flow

        self.robust_pressure_load              = tonp(tM.addVars(self.gas_load_num, self.T))                        # this is unknown load/gen node pressure
        self.robust_pressure_generator         = tonp(tM.addVars(self.gas_generator_num, self.T))                        # this is unknown load/gen node pressure

        self.zUp                               = tonp(tM.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zDown                             = tonp(tM.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zUp[0,0] = 1
        self.zUp[0,1] = 0
        self.zDown[0,0] = 0
        self.zDown[0,1] = 1

        # build and test_pg feasible model
        self.feasibleProblemNetworkBuildTest(self.testModel)
    def feasibleProblemNetworkBuildTest(self, robustModel):
        dG = self.dualGen

        # build network constraints
        MldPsr = np.vstack((
            self.robust_flow_load                 .reshape(-1, 1),         # this is known value
            self.robust_gas_generator_trans       .reshape(-1, 1),         # this is var
            self.robust_pressure_well             .reshape(-1, 1)          # this is known value
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
        #
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

        # generate network constraints
        print('----generate robust Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count - byPassIndex):
            row = row + byPassIndex
            if row % 100 == 0:
                print('----generate robust Step 1 ' + str(row) + ' ----------------')
                print('common')
            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     sense='=='
            # )
            # continue
            coeffInput = [ -1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (1 * np.sum(coeff, axis=1)).tolist())

            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex. extend( MAGIC_LIST )

            constantInput = 1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '<=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

            coeffInput = [ 1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (-1 * np.sum(coeff, axis=1)).tolist())
            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex.extend( MAGIC_LIST )

            constantInput = -1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '>=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

        print('----generate robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                dG.addConstr(
                    self.robust_gas_generator[gen, t],
                    (self.robust_gas_generator_base[gen, t] +
                    self.zUp[gen, t] * self.robust_gas_generator_reserve_up[gen, t] -
                    self.zDown[gen, t] * self.robust_gas_generator_reserve_down[gen, t]),
                    sense='==',
                    appendIndex=True
                )

        print('----generate robust Step 3----------------')
        # build pressure limit
        sCollection = []
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

        # for gas_gen in self.robust_gas_generator.flatten().tolist():
        #     dG.addConstr(
        #         gas_gen, 0, sense='>='
        #     )


        print('----generate robust Step 3.1 add objective----------------')
        # build primitive objective
        dG.addObjectiveMin(gurobi.quicksum(sCollection))

        print('----generate robust Step 3.2 get dual problem----------------')
        # add additional objective
        dualModel = dG.getDual()

        # print('----generate robust Step 3.3 add additional objective----------------')
        # dualIndexList = dG.robustIndex
        # dualIndexCount = 0
        # BIGM = 1e3
        # oobbjj = []
        #
        #
        # self.zUp = tonp(dualModel.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        # self.zDown = tonp(dualModel.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        # self.zGreatV1 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        # self.zGreatV2 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        # self.zLessV1 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        # self.zLessV2 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        #
        #
        # for gen in range(self.gas_generator_num):
        #     for t in range(self.T_long):    # per equal constraints
        #         dual1 = dG.dualVars[dualIndexList[dualIndexCount]]
        #         dualIndexCount = dualIndexCount + 1
        #         dual2 = dG.dualVars[dualIndexList[dualIndexCount]]
        #         dualIndexCount = dualIndexCount + 1
        #
        #         assert not dual1.sameAs(dual2)
        #
        #         zUp = self.zUp[gen, int(t/self.time_per_T)]
        #         zDown = self.zDown[gen, int(t/self.time_per_T)]
        #         dualModel.addConstr(zUp + zDown <= 1)
        #
        #         # for great part
        #         gasBase = self.robust_gas_generator_base[gen, t]
        #         gasUpper = self.robust_gas_generator_reserve_up[gen, t]
        #         gasDown = self.robust_gas_generator_reserve_down[gen, t]
        #         # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
        #         v1 = dualModel.addVar(lb=-1*INF, ub=0)
        #         v2 = dualModel.addVar(lb=-1*INF, ub=0)
        #         self.zGreatV1[gen, t] = v1
        #         self.zGreatV2[gen, t] = v2
        #         # linear     v1 = zUp * dual1
        #         self.addRelaxConstraints(v1, zUp, dual1, BIGM, dualModel)
        #         # linear    v2 = zDown * dual1
        #         self.addRelaxConstraints(v2, zDown, dual1, BIGM, dualModel)
        #         oobbjj.append(dual1 * gasBase + gasUpper * v1 - gasDown * v2)
        #
        #         # for less part
        #         # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
        #         v1 = dualModel.addVar(lb=-1*INF, ub=0)
        #         v2 = dualModel.addVar(lb=-1*INF, ub=0)
        #         self.zLessV1[gen, t] = v1
        #         self.zLessV2[gen, t] = v2
        #         # linear     v1 = zUp * dual2
        #         self.addRelaxConstraints(v1, zUp, dual2, BIGM, dualModel)
        #         # linear    v2 = zDown * dual2
        #         self.addRelaxConstraints(v2, zDown, dual2, BIGM, dualModel)
        #         oobbjj.append(-1 * dual2 * gasBase - gasUpper * v1 + gasDown * v2)
        #
        # assert dualIndexCount == len(dualIndexList)


        print('----generate robust Step 3.4 optimize----------------')
        # dualModel.setObjective( 1*dualModel.getObjective() + 0*sum(oobbjj) )

        # dualModel.setParam("MIPFocus", 1)
        dualModel.optimize()

        self.zGreat = to_value(self.zUp)
        self.zLess = to_value(self.zDown)

        return dualModel.getObjective().getValue()





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
    def buildRobustVarsNetwork(self, model):
        self.robust_pressure_well              = to_value(self.node_pressure_trans[0])                         # this is known source pressure

        self.robust_gas_generator_base         = to_value(self.gas_generator)                                                          # this is known generator load
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
    def addRelaxConstraints(self, v, z, dual, BIGM, model):
        model.addConstr(-1 * BIGM * z <= v)
        model.addConstr(v <= 0)
        model.addConstr(-1 * BIGM * (1 - z) <= (dual - v))
        model.addConstr((dual - v) <= 0)
    def feasibleProblemNetworkPreset(self):
        # build robust var
        self.robustModel = gurobi.Model()
        self.dualGen = GenDual(self.robustModel)

        rM = self.robustModel
        self.buildRobustVarsNetwork(rM)
    def feasibleProblemNetworkBuild(self, robustModel):
        dG = self.dualGen

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
            # continue

            coeffInput = [ -1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_fre, -1))
            coeffInput.extend( (1 * np.sum(coeff, axis=1)).tolist())

            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex. extend( MAGIC_LIST )

            constantInput = 1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            assert(len(coeffInput) == len(varIndex))
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '<=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

            coeffInput = [ 1]
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_fre, -1))
            coeffInput.extend( (-1 * np.sum(coeff, axis=1)).tolist())
            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex.extend( MAGIC_LIST )

            constantInput = -1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '>=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

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
        sCollection = []
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

        assert dualIndexCount == len(dualIndexList)
        dualModel.setObjective( 1*dualModel.getObjective() + 1*sum(oobbjj) )

        if LOG_OUTPUT:
            print('----generate robust Step 3.4 optimize----------------')
        # dualModel.setParam("MIPFocus", 1)
        dualModel.setParam("OutputFlag", 0)
        dualModel.optimize()

        self.zGreat = to_value(self.zUp)
        self.zLess = to_value(self.zDown)

        return dualModel.getObjective().getValue()






    # [DEPRECATED] this is original feasible problem with network form
    # but we change to function to two sub-function as :
    # preset and build as shown above, the decouple the
    # var and function build.
    def feasibleProblemNetwork(self):
        # build robust var
        self.robustModel = gurobi.Model()
        self.dualGen = GenDual(self.robustModel)

        rM = self.robustModel
        dG = self.dualGen

        self.buildRobustVarsNetwork(rM)
        # build network constraints
        MldPsr = np.vstack((
            self.robust_flow_load           .reshape(-1, 1),         # this is known value
            self.robust_gas_generator_trans       .reshape(-1, 1),         # this is var
            self.robust_pressure_well       .reshape(-1, 1)          # this is known value
        ))
        start_index = self.robust_flow_load.reshape(-1, 1).shape[0]
        middle_index = self.robust_flow_load.reshape(-1, 1).shape[0] + self.robust_gas_generator_trans.reshape(-1, 1).shape[0]

        MsrPld = np.vstack((
            self.robust_flow_source        .reshape(-1, 1),
            self.robust_pressure_load      .reshape(-1, 1),
            self.robust_pressure_generator .reshape(-1, 1)       # all these are vars
        ))

        # because of the O^2 complex of index search , so we pre-set variable index
        dG.preset_var_index((
            self.robust_gas_generator.reshape(-1, 1).flatten().tolist(),
            self.robust_flow_source.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist()
        ))
        MAGIC_LIST = np.arange(self.robust_gas_generator.reshape(-1, 1).shape[0]).tolist()

        map_value_index = 0
        index_m = {}
        self.robustModel.update()
        for vars_preset in (
            self.robust_gas_generator.reshape(-1, 1).flatten().tolist(),
            self.robust_flow_source.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_load.reshape(-1, 1).flatten().tolist(),
            self.robust_pressure_generator.reshape(-1, 1).flatten().tolist()
        ):
            vars_list = vars_preset
            for var in vars_list:
                index_m[var] = map_value_index
                map_value_index = map_value_index + 1

        # MAGIC_LIST = [index_m[gas_gen[0]]
        #               for gas_gen in self.robust_gas_generator_trans.reshape(-1, 1).tolist()]

        print('----generate robust Step 1 ' + str(self.YY.shape[0]) + '----------------')
        iter_count = self.YY.shape[0]
        for row in range(iter_count):
            if row % 100 == 0:
                print('----generate robust Step 1 ' + str(row) + ' ----------------')
            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     sense='=='
            # )

            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     sense='<='
            # )
            # MsrPld[row, 0] - (self.YY[row, :]).dot(MldPsr)[0] <= 0   aka C
            #                   self.YY[row, :] <====> [ .... start ..... middle ..... ]
            # the value part is :
            # self.YY[row, :start]. dot( MldPsr[:start] ) +
            # self.YY[row, middle:].dot( MldPsr[middle:] )      aka A
            # the var part is :
            # self.YY[start:middle].dot( MldPsr[start:middle] )  aka B
            # so the C is equal as:
            # MsrPld[row, 0] - ( A + B )   ===> MsrPld[row, 0] - B - A
            # =====> [ 1, -self.YY[row, start:middle] ] * [ MsrPld[row, 0], MldPsr[start:middle] ] - A <= 0
            # =====> [ -1, self.YY[row, start:middle] ] * [ MsrPld[row, 0], MldPst[start:middle] ] + A >= 0
            # so:
            # coeffInput=None, varInput=None, constantInput=None
            coeffInput = [ -1]
            # coeffInput.extend( (1 * np.array(self.YY[row, start_index:middle_index])).tolist() )
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (1 * np.sum(coeff, axis=1)).tolist())

            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex. extend( MAGIC_LIST )

            constantInput = 1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '<=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

            # dG.addConstr(
            #     MsrPld[row, 0],
            #     (self.YY[row, :]).dot(MldPsr)[0],
            #     sense='>='
            # )
            # MsrPld[row, 0] - (self.YY[row, :]).dot(MldPsr)[0] >= 0   aka C
            #                   self.YY[row, :] <====> [ .... start ..... middle ..... ]
            # the value part is :
            # self.YY[row, :start]. dot( MldPsr[:start] ) +
            # self.YY[row, middle:].dot( MldPsr[middle:] )      aka A
            # the var part is :
            # self.YY[row, start:middle].dot( MldPsr[start:middle] )  aka B
            # so the C is equal as:
            # MsrPld[row, 0] - ( A + B )   ===> MsrPld[row, 0] - B - A
            # =====> [ 1, -self.YY[row, start:middle] ] * [ MsrPld[row, 0], MldPsr[start:middle] ] - A >= 0
            # so:
            # coeffInput=None, varInput=None, constantInput=None
            coeffInput = [ 1]
            # coeffInput.extend( (-1 * np.array(self.YY[row, start_index:middle_index])).tolist() )
            coeff = self.YY[row, start_index:middle_index]
            coeff = coeff.reshape((self.T_long, -1))
            coeffInput.extend( (-1 * np.sum(coeff, axis=1)).tolist())
            varInput =   [ MsrPld[row, 0] ]
            varInput.extend( self.robust_gas_generator.flatten().tolist() )

            varIndex =   [ index_m[MsrPld[row, 0]]]
            varIndex.extend( MAGIC_LIST )

            constantInput = -1 * (
                self.YY[row, :start_index].dot(MldPsr[:start_index])[0] +
                self.YY[row, middle_index:].dot(MldPsr[middle_index:])[0]
            )
            dG.addConstr(
                MsrPld[row, 0],
                (self.YY[row, :]).dot(MldPsr)[0],
                '>=',
                name='xx',
                appendIndex=False,
                coeffInput=coeffInput,
                varInput=varInput,
                constantInput=constantInput,
                varIndex=varIndex,
                showTime=False
            )

        print('----generate robust Step 2----------------')
        # build gas generator output
        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):
                dG.addConstr(
                    self.robust_gas_generator[gen, t],
                    0,
                    # (self.robust_gas_generator_base[gen, t] +
                    # self.zUp[gen, int(t/self.time_per_T)] * self.robust_gas_generator_reserve_up[gen, int(t/self.time_per_T)] -
                    # self.zDown[gen, int(t/self.time_per_T)] * self.robust_gas_generator_reserve_down[gen, int(t/self.time_per_T)]),
                    sense='==',
                    appendIndex=True
                )

        print('----generate robust Step 3----------------')
        # build pressure limit
        sCollection = []
        node_collection = self.robust_pressure_load.flatten().tolist() + \
                          self.robust_pressure_generator.flatten().tolist()
        for load in node_collection:
            sMin = self.robustModel.addVar()
            sPlus = self.robustModel.addVar()
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


        # build primitive objective
        dG.addObjectiveMin(gurobi.quicksum(sCollection))

        # add additional objective
        dualModel = dG.getDual()

        dualIndexList = dG.robustIndex
        dualIndexCount = 0
        BIGM = 1e3
        oobbjj = []


        self.zUp = tonp(dualModel.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zDown = tonp(dualModel.addVars(self.gas_generator_num, self.T_long, vtype=gurobi.GRB.BINARY))
        self.zGreatV1 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        self.zGreatV2 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        self.zLessV1 = np.empty((self.gas_generator_num, self.T), dtype=np.object)
        self.zLessV2 = np.empty((self.gas_generator_num, self.T), dtype=np.object)


        for gen in range(self.gas_generator_num):
            for t in range(self.T_long):    # per equal constraints
                dual1 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1
                dual2 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1

                assert not dual1.sameAs(dual2)

                zUp = self.zUp[gen, int(t/self.time_per_T)]
                zDown = self.zDown[gen, int(t/self.time_per_T)]
                dualModel.addConstr(zUp + zDown <= 1)

                # for great part
                gasBase = self.robust_gas_generator_base[gen, t]
                gasUpper = self.robust_gas_generator_reserve_up[gen, t]
                gasDown = self.robust_gas_generator_reserve_down[gen, t]
                # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
                v1 = dualModel.addVar(lb=-1*INF, ub=0)
                v2 = dualModel.addVar(lb=-1*INF, ub=0)
                self.zGreatV1[gen, t] = v1
                self.zGreatV2[gen, t] = v2
                # v1 = zUp * dual
                # v2 = zDown * dual
                # temp1 = gasRobust * dual = gasBase * dual + zUp * gasUpper * dual - zDown * gasDown * dual
                # temp2                    = gasBase * dual + gasUpper * v1 - gasDown * v2

                # dualModel.addConstr(-1 * BIGM * zUp <= v1,                  )
                # dualModel.addConstr(v1 <= 0,                                      )
                # dualModel.addConstr(-1 * BIGM * (1 - zUp) <= (dual1 - v1),  )
                # dualModel.addConstr((dual1 - v1) <= 0,                            )

                # linear     v1 = zUp * dual1
                self.addRelaxConstraints(v1, zUp, dual1, BIGM, dualModel)
                # dualModel.addConstr(-1 * BIGM * zDown <= v2,                  )
                # dualModel.addConstr(v2 <= 0,                                      )
                # dualModel.addConstr(-1 * BIGM * (1 - zDown) <= (dual1 - v2),  )
                # dualModel.addConstr((dual1 - v2) <= 0,                            )

                # linear    v2 = zDown * dual1
                self.addRelaxConstraints(v2, zDown, dual1, BIGM, dualModel)
                # dualModel.addConstr(v1 <= 0)
                # dualModel.addConstr(dual <= v1)
                # dualModel.addConstr(-1 * z1 <= v1)
                # dualModel.addConstr(v1 <= -1 * z1 + dual + 1)
                #
                # dualModel.addConstr(v2 <= 0)
                # dualModel.addConstr(dual <= v2)
                # dualModel.addConstr(-1 * z2 <= v2)
                # dualModel.addConstr(v2 <= -1 * z2 + dual + 1)
                oobbjj.append(dual1 * gasBase + gasUpper * v1 - gasDown * v2)

                # for less part
                # gasRobust = gasBase + zUp * gasUpper - zDown * gasDown
                v1 = dualModel.addVar(lb=-1*INF, ub=0)
                v2 = dualModel.addVar(lb=-1*INF, ub=0)
                self.zLessV1[gen, t] = v1
                self.zLessV2[gen, t] = v2
                # v1 = zUp * dual
                # v2 = zDown * dual
                # temp1 = gasRobust * dual = gasBase * dual + zUp * gasUpper * dual - zDown * gasDown * dual
                # temp2                    = gasBase * dual + gasUpper * v1 - gasDown * v2

                # linear     v1 = zUp * dual2
                # dualModel.addConstr(-1 * BIGM * zUp <= v1                    )
                # dualModel.addConstr(v1 <= 0                                  )
                # dualModel.addConstr(-1 * BIGM * (1 - zDown) <= (dual2 - v1)  )
                # dualModel.addConstr((dual2 - v1) <= 0                        )

                # linear     v1 = zUp * dual2
                self.addRelaxConstraints(v1, zUp, dual2, BIGM, dualModel)
                # linear    v2 = zDown * dual2
                # dualModel.addConstr(-1 * BIGM * zDown <= v2                 )
                # dualModel.addConstr(v2 <= 0                                 )
                # dualModel.addConstr(-1 * BIGM * (1 - zDown) <= (dual2 - v2) )
                # dualModel.addConstr((dual2 - v2) <= 0                       )

                # linear    v2 = zDown * dual2
                self.addRelaxConstraints(v2, zDown, dual2, BIGM, dualModel)
                # dualModel.addConstr(v1 <= 0)
                # dualModel.addConstr(dual <= v1)
                # dualModel.addConstr(-1 * z1 <= v1)
                # dualModel.addConstr(v1 <= -1 * z1 + dual + 1)
                #
                # dualModel.addConstr(v2 <= 0)
                # dualModel.addConstr(dual <= v2)
                # dualModel.addConstr(-1 * z2 <= v2)
                # dualModel.addConstr(v2 <= -1 * z2 + dual + 1)
                oobbjj.append(-1 * dual2 * gasBase - gasUpper * v1 + gasDown * v2)

        # MAGIC_NUM = int(10)
        # for node in range(self.gas_node_num):
        #     for gas in range(self.gas_generator_num):
        #         for t in range(self.T):
        #             dualModel.addConstr(
        #                 self.zGreat[node, t, gas] ==
        #                 self.zGreat[node, int(t/MAGIC_NUM) * MAGIC_NUM, gas])
        #             dualModel.addConstr(
        #                 self.zLess[node, t, gas] ==
        #                 self.zLess[node, int(t/MAGIC_NUM) * MAGIC_NUM, gas])

        assert dualIndexCount == len(dualIndexList)
        dualModel.setObjective( dualModel.getObjective() + sum(oobbjj) )

        # dualModel.setParam('OutputFlag', 0)
        # # print("=====> statistic start ")
        # dualModel.update()
        # dualModel.write(f'RobustModel_{self.file_index}.lp')
        # self.file_index = self.file_index + 1
        # # print("=====> statistic end ")
        # dualModel.tune()
        # dualModel.setParam("MIPFocus", 1)
        dualModel.optimize()

        self.zGreat = to_value(self.zUp)
        self.zLess = to_value(self.zDown)

        return dualModel.getObjective().getValue()
    def buildRobustVars(self, robustModel):
        self.robust_start_pipeline_flow_trans = tonp(robustModel.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_flow_min, ub=self.gas_flow_max,
            name='pipeline_gas_flow_of_start_port'
        ))
        self.robust_end_pipeline_flow_trans = tonp(robustModel.addVars(
            self.gas_line_num, self.T,
            lb=self.gas_flow_min, ub=self.gas_flow_max,
            name='pipeline_gas_flow_of_end_port'
        ))
        self.robust_node_pressure_trans = tonp(robustModel.addVars(
            self.gas_node_num, self.T,
            lb=self.gas_node_pressure_min, ub=self.gas_node_pressure_max,
            name='gas_node_pressure_trans'
        ))
        self.robust_gas_generator_trans = tonp(robustModel.addVars(
            self.gas_generator_num, self.T,
            name='gas_generator_output_trans'
        ))

        self.robust_well_in_trans = tonp(robustModel.addVars(
            self.gas_well_num, self.T,
            name='gas_well_in_trans_robust'
        ))

        self.robust_gas_load_trans = to_value(self.gas_load_trans)
        # self.robust_well_in_trans = to_value(self.well_in_trans)
        self.first_stage_node_pressure_trans = to_value(self.node_pressure_trans)
    # [DEPRECATED] this is feasible problem without network form
    def feasibleProblem(self):
        self.robustModel = gurobi.Model()
        self.dualGen = GenDual(self.robustModel)

        dG = self.dualGen
        rM = self.robustModel

        self.buildRobustVars(rM)

        def myDot2(matrix, syms, row, y_ranges):
            value = matrix[row][self.node_counts - 1][y_ranges]
            sym = syms[y_ranges]
            return sum(value * sym)

        # print('======> Robust finish 7.1')
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

                start_pipeline_flow = self.robust_start_pipeline_flow_trans[line]
                end_pipeline_flow = self.robust_end_pipeline_flow_trans[line]
                start_pipeline_pressure = self.robust_node_pressure_trans[start_node]
                end_pipeline_pressure = self.robust_node_pressure_trans[end_node]

                flow2pressure1 = np.hstack((end_pipeline_flow, start_pipeline_pressure))

                # print('=======> finish 7.1.1')
                dG.addConstr(start_pipeline_flow[t],
                             myDot2(self.Matrix_ab, flow2pressure1, t, y_range),
                             sense='==',
                             name='ab[' + str(line) + ',' + str(t) + ']'
                             )
                dG.addConstr(end_pipeline_pressure[t],
                             myDot2(self.Matrix_cd, flow2pressure1, t, y_range2),
                             sense='==',
                             name='cd[' + str(line) + ',' + str(t) + ']'
                             )

        # print('======> Robust finish 7.2')
        self.sCollection = []
        for t in range(self.T):
            for node in range(self.gas_node_num):
                sPlus = rM.addVar()
                sMinus = rM.addVar()
                self.sCollection.append(sPlus)
                self.sCollection.append(sMinus)

                dG.addConstr(
                    exprLeft=
                    sum(self.robust_end_pipeline_flow_trans[  np.where(self.gas_pipeline_end_node == node),     t].flatten()) +
                    sum(self.robust_well_in_trans[            np.where(self.gas_well_connection_index == node), t].flatten()) -
                    sum(self.robust_start_pipeline_flow_trans[np.where(self.gas_pipeline_start_node == node),   t].flatten()) +
                    sPlus - sMinus,
                    exprRight=
                    sum(self.robust_gas_load_trans[           np.where(self.gas_load_connection_index == node), t].flatten()),
                    sense='==',
                    appendIndex=True,
                    showTime=False,
                    name='nodeBalance[' + str(node) + ',' + str(t) + ']'
                )
                dG.addConstr(sPlus, 0, '>=', name='sPZero[' + str(node) + ',' + str(t) + ']')
                dG.addConstr(sMinus, 0, '>=', name='sMZero[' + str(node) + ',' + str(t) + ']')

        # ====> node pressure specify [gas system transient]
        # print('======> Robust finish 7.3')
        for t in range(self.T):
            for well in range(self.gas_well_num):
                node = self.gas_well_connection_index[well]
                dG.addConstr(
                    self.robust_node_pressure_trans[node, t],
                    self.first_stage_node_pressure_trans[node, t],
                    sense='==',
                    name='WellPress[' + str(well) + ',' + str(t) + ']'
                )

        # =====> node pressure max/min limit
        # print('======> Robust finish 7.4')
        for t in range(self.T):
            for node in range(self.gas_node_num):
                dG.addConstr(
                    self.robust_node_pressure_trans[node, t],
                    self.gas_node_pressure_min,
                    sense='>=',
                    name='nodePressMin[' + str(node) + ',' + str(t) + ']'
                )

                dG.addConstr(
                    self.robust_node_pressure_trans[node, t],
                    self.gas_node_pressure_max,
                    sense='<=',
                    name='nodePressMax[' + str(node) + ',' + str(t) + ']'
                )

        for t in range(self.T):
            for line in range(self.gas_line_num):
                dG.addConstr(
                    self.robust_start_pipeline_flow_trans[line, t],
                    self.gas_flow_min,
                    sense='>=',
                    name='startFlowMin[' + str(line) + ',' + str(t) + ']'
                )

                dG.addConstr(
                    self.robust_start_pipeline_flow_trans[line, t],
                    self.gas_flow_max,
                    sense='<=',
                    name='startFlowMax[' + str(line) + ',' + str(t) + ']'
                )

                dG.addConstr(
                    self.robust_end_pipeline_flow_trans[line, t],
                    self.gas_flow_min,
                    sense='>=',
                    name='endFlowMin[' + str(line) + ',' + str(t) + ']'
                )

                dG.addConstr(
                    self.robust_end_pipeline_flow_trans[line, t],
                    self.gas_flow_max,
                    sense='<=',
                    name='endFlowMax[' + str(line) + ',' + str(t) + ']'
                )

        # print('======> Robust finish 7.5')
        dG.addObjectiveMin(sum(self.sCollection) * 1)
        dualModel = dG.getDual()


        # print('======> Robust finish 7.6')
        # add additional robust objective part
        dualIndexList = dG.robustIndex
        dualIndexCount = 0
        BIGM = 1e3
        oobbjj = []

        self.zGreat = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)
        self.zLess = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)
        self.zGreatV1 = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)
        self.zGreatV2 = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)
        self.zLessV1 = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)
        self.zLessV2 = np.empty((self.gas_node_num, self.T, self.gas_generator_num), dtype=np.object)


        for t in range(self.T):
            for node in range(self.gas_node_num):    # per equal constraints
                dual1 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1
                dual2 = dG.dualVars[dualIndexList[dualIndexCount]]
                dualIndexCount = dualIndexCount + 1

                assert dualIndexList[dualIndexCount - 2] + 1 == dualIndexList[dualIndexCount - 1]
                assert not dual1.sameAs(dual2)

                # 这是当前时刻，连接该节点的所有燃气轮机的 出力 的 集合
                expr1 = self.robust_gas_generator_trans[np.where(self.gas_gen_index_gas == node), t].flatten() * 1
                # 这是但钱时刻，连接该节点的所有燃气轮机的 出力 的负值 的集合
                expr2 = self.robust_gas_generator_trans[np.where(self.gas_gen_index_gas == node), t].flatten() * -1
                # 这是 第几个 燃气轮机连接在 当前节点的 集合
                windIndexList = np.where(self.gas_gen_index_gas == node)[0].tolist()

                z1 = dualModel.addVars(len(expr1), vtype=gurobi.GRB.BINARY, name='gasUpBin')
                z2 = dualModel.addVars(len(expr1), vtype=gurobi.GRB.BINARY, name='gasDownBin')
                for i in range(len(expr1)):
                    dualModel.addConstr(z1[i] + z2[i] <= 1)

                # for great part
                for index, item in enumerate(expr1):
                    self.zGreat[node, t, index] = z1[index]

                    gasIndex = windIndexList[index]
                    gasBase = self.gas_generator_trans[gasIndex, int(t/self.time_per_T)].getAttr('X')
                    gasUpper = self.gas_generator_reserve_up[gasIndex, int(t/self.time_per_T)].getAttr('X')
                    gasDown = self.gas_generator_reserve_down[gasIndex, int(t/self.time_per_T)].getAttr('X')

                    # z1 = dualModel.addVar(type=gurobi.GRB.BINARY)
                    # z2 = dualModel.addVar(type=gurobi.GRB.BINARY)

                    # gasRobust = gasBase + z1[index] * gasUpper - z2[index] * gasDown
                    # dualModel.addConstr(z1 + z2 <= 1)

                    v1 = dualModel.addVar(lb=-1*INF, ub=0, name='gasUpAuxGreat['   + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    v2 = dualModel.addVar(lb=-1*INF, ub=0, name='gasDownAuxGreat[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    self.zGreatV1[node, t, index] = v1
                    self.zGreatV2[node, t, index] = v2
                    # v1 = z1[index] * dual
                    # v2 = z2[index] * dual
                    # temp1 = gasRobust * dual = gasBase * dual + z1[index] * gasUpper * dual - z2[index] * gasDown * dual
                    # temp2                    = gasBase * dual + gasUpper * v1 - gasDown * v2

                    # linear     v1 = z1 * dual1
                    dualModel.addConstr(-1 * BIGM * z1[index] <= v1,                  name='relaxGreat1[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(v1 <= 0,                                      name='relaxGreat2[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(-1 * BIGM * (1 - z1[index]) <= (dual1 - v1),  name='relaxGreat3[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr((dual1 - v1) <= 0,                            name='relaxGreat4[' + str(node) + ',' + str(index) + ',' + str(t) + ']')

                    # linear    v2 = z2 * dual1
                    dualModel.addConstr(-1 * BIGM * z2[index] <= v2,                  name='relaxGreat5[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(v2 <= 0,                                      name='relaxGreat6[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(-1 * BIGM * (1 - z2[index]) <= (dual1 - v2),  name='relaxGreat7[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr((dual1 - v2) <= 0,                            name='relaxGreat8[' + str(node) + ',' + str(index) + ',' + str(t) + ']')

                    # dualModel.addConstr(v1 <= 0)
                    # dualModel.addConstr(dual <= v1)
                    # dualModel.addConstr(-1 * z1 <= v1)
                    # dualModel.addConstr(v1 <= -1 * z1 + dual + 1)
                    #
                    # dualModel.addConstr(v2 <= 0)
                    # dualModel.addConstr(dual <= v2)
                    # dualModel.addConstr(-1 * z2 <= v2)
                    # dualModel.addConstr(v2 <= -1 * z2 + dual + 1)
                    oobbjj.append(dual1 * gasBase + gasUpper * v1 - gasDown * v2)

                # for less part
                for index, item in enumerate(expr1):
                    self.zLess[node, t, index] = z2[index]

                    gasIndex = windIndexList[index]
                    gasBase = self.gas_generator_trans[gasIndex, int(t/self.time_per_T)].getAttr('X')
                    gasUpper = self.gas_generator_reserve_up[gasIndex, int(t/self.time_per_T)].getAttr('X')
                    gasDown = self.gas_generator_reserve_down[gasIndex, int(t/self.time_per_T)].getAttr('X')

                    # z1 = dualModel.addVar(type=gurobi.GRB.BINARY)
                    # z2 = dualModel.addVar(type=gurobi.GRB.BINARY)

                    # gasRobust = gasBase + z1[index] * gasUpper - z2[index] * gasDown
                    # dualModel.addConstr(z1 + z2 <= 1)

                    v1 = dualModel.addVar(lb=-1*INF, ub=0, name='gasUpAuxLess['   + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    v2 = dualModel.addVar(lb=-1*INF, ub=0, name='gasDownAuxLess[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    self.zLessV1[node, t, index] = v1
                    self.zLessV2[node, t, index] = v2
                    # v1 = z1[index] * dual
                    # v2 = z2[index] * dual
                    # temp1 = gasRobust * dual = gasBase * dual + z1[index] * gasUpper * dual - z2[index] * gasDown * dual
                    # temp2                    = gasBase * dual + gasUpper * v1 - gasDown * v2

                    # linear     v1 = z1 * dual2
                    dualModel.addConstr(-1 * BIGM * z1[index] <= v1,                 name='relaxLess1[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(v1 <= 0,                                     name='relaxLess2[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(-1 * BIGM * (1 - z1[index]) <= (dual2 - v1), name='relaxLess3[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr((dual2 - v1) <= 0,                           name='relaxLess4[' + str(node) + ',' + str(index) + ',' + str(t) + ']')

                    # linear    v2 = z2 * dual2
                    dualModel.addConstr(-1 * BIGM * z2[index] <= v2,                 name='relaxLess5[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(v2 <= 0,                                     name='relaxLess6[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr(-1 * BIGM * (1 - z2[index]) <= (dual2 - v2), name='relaxLess7[' + str(node) + ',' + str(index) + ',' + str(t) + ']')
                    dualModel.addConstr((dual2 - v2) <= 0,                           name='relaxLess8[' + str(node) + ',' + str(index) + ',' + str(t) + ']')

                    # dualModel.addConstr(v1 <= 0)
                    # dualModel.addConstr(dual <= v1)
                    # dualModel.addConstr(-1 * z1 <= v1)
                    # dualModel.addConstr(v1 <= -1 * z1 + dual + 1)
                    #
                    # dualModel.addConstr(v2 <= 0)
                    # dualModel.addConstr(dual <= v2)
                    # dualModel.addConstr(-1 * z2 <= v2)
                    # dualModel.addConstr(v2 <= -1 * z2 + dual + 1)
                    oobbjj.append(-1 * dual2 * gasBase - gasUpper * v1 + gasDown * v2)

        MAGIC_NUM = int(10)
        for node in range(self.gas_node_num):
            for gas in range(self.gas_generator_num):
                for t in range(self.T):
                    dualModel.addConstr(
                        self.zGreat[node, t, gas] ==
                        self.zGreat[node, int(t/MAGIC_NUM) * MAGIC_NUM, gas])
                    dualModel.addConstr(
                        self.zLess[node, t, gas] ==
                        self.zLess[node, int(t/MAGIC_NUM) * MAGIC_NUM, gas])

        assert dualIndexCount == len(dualIndexList)
        dualModel.setObjective( dualModel.getObjective() + sum(oobbjj) )

        # dualModel.setParam('OutputFlag', 0)
        # # print("=====> statistic start ")
        # dualModel.update()
        # dualModel.write(f'RobustModel_{self.file_index}.lp')
        # self.file_index = self.file_index + 1
        # # print("=====> statistic end ")
        # dualModel.tune()
        # dualModel.setParam("MIPFocus", 1)
        dualModel.optimize()

        self.zGreat = to_value(self.zGreat)
        self.zLess = to_value(self.zLess)

        return dualModel.getObjective().getValue()







    #
    def addAppendVars(self):
        self.append_pressure_well         = self.node_pressure_trans[0, :]
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
    # [DEPRECATED]this is append constraint without network form
    def appendConstraint(self):
        # self.robustModel = gurobi.Model()
        # self.dualGen = GenDual(self.robustModel)
        #
        # dG = self.dualGen
        # rM = self.robustModel

        self.buildRobustVars(self.model)

        def myDot2(matrix, syms, row, y_ranges):
            value = matrix[row][self.node_counts - 1][y_ranges]
            sym = syms[y_ranges]
            return sum(value * sym)

        # print('======> Robust finish 7.1')
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

                start_pipeline_flow = self.robust_start_pipeline_flow_trans[line]
                end_pipeline_flow = self.robust_end_pipeline_flow_trans[line]
                start_pipeline_pressure = self.robust_node_pressure_trans[start_node]
                end_pipeline_pressure = self.robust_node_pressure_trans[end_node]

                flow2pressure1 = np.hstack((end_pipeline_flow, start_pipeline_pressure))

                # print('=======> finish 7.1.1')
                self.model.addConstr(start_pipeline_flow[t] ==
                                     myDot2(self.Matrix_ab, flow2pressure1, t, y_range)
                                     )
                self.model.addConstr(end_pipeline_pressure[t] ==
                                     myDot2(self.Matrix_cd, flow2pressure1, t, y_range2)
                                     )

        # print('======> Robust finish 7.2')
        for t in range(self.T):
            for node in range(self.gas_node_num):
                self.model.addConstr(
                    sum(self.robust_end_pipeline_flow_trans[  np.where(self.gas_pipeline_end_node == node),     t].flatten()) +
                    sum(self.robust_well_in_trans[            np.where(self.gas_well_connection_index == node), t].flatten()) -
                    sum(self.robust_start_pipeline_flow_trans[np.where(self.gas_pipeline_start_node == node),   t].flatten())
                    ==
                    sum(self.gas_load_trans[                  np.where(self.gas_load_connection_index == node), t].flatten()) +
                    sum((self.gas_generator_trans[np.where(self.gas_gen_index_gas == node), t] +
                        self.zGreat[node, t] * self.gas_generator_reserve_up[np.where(self.gas_gen_index_gas == node), int(t/self.time_per_T)] -
                        self.zLess[node, t ] * self.gas_generator_reserve_down[np.where(self.gas_gen_index_gas == node), int(t/self.time_per_T)]).flatten())
                )

        # ====> node pressure specify [gas system transient]
        # print('======> Robust finish 7.3')
        for t in range(self.T):
            for well in range(self.gas_well_num):
                node = self.gas_well_connection_index[well]

                self.model.addConstr(
                    self.robust_node_pressure_trans[node, t] ==
                    self.node_pressure_trans[node, t]
                )

        # =====> node pressure max/min limit
        # print('======> Robust finish 7.4')
        for t in range(self.T):
            for node in range(self.gas_node_num):
                self.model.addConstr(
                    self.robust_node_pressure_trans[node, t] >=
                    self.gas_node_pressure_min,
                )

                self.model.addConstr(
                    self.robust_node_pressure_trans[node, t] <=
                    self.gas_node_pressure_max,
                )

        for t in range(self.T):
            for line in range(self.gas_line_num):
                self.model.addConstr(
                    self.robust_start_pipeline_flow_trans[line, t] >=
                    self.gas_flow_min,
                )

                self.model.addConstr(
                    self.robust_start_pipeline_flow_trans[line, t] <=
                    self.gas_flow_max,
                )

                self.model.addConstr(
                    self.robust_end_pipeline_flow_trans[line, t] >=
                    self.gas_flow_min,
                )

                self.model.addConstr(
                    self.robust_end_pipeline_flow_trans[line, t] <=
                    self.gas_flow_max,
                 )

        self.model.optimize()












    def drawPressure(self):
        import matplotlib.pyplot as plt

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
    pg.buildBasePrimaryProblem()

    # pg.model.optimize()
    # pg.feasibleProblem()
    # pg.appendConstraint()


    for i in range(17):
        print('Stage 1 : preSolve : [ Iteration: ' + str(i) + ']')
        pg.model.setParam('OutputFlag', 0)
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
        t2 = time.time()
        if LOG_OUTPUT:
            print('Feasibel Time : ' + str(t2 - t1))
        if abs(feasibleTest) <= 1e-1:
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
            pg.appendConstraintNetwork()
        print(' ')



def doTest():
    # 1, test_pg network form generation correction
    from test_pg.network6node import getConfig as testConfig
    pg_test = PowerGas(*testConfig())



if __name__ == '__main__':
    # testGetDual()
    # assert 0

    pg = PowerGas(*getConfig())

    # pg.testMatrix()
    # assert 0


    scenarioLess = []
    scenarioGreat = []

    doMain()
    # doTest()


def drawScenario():
    import matplotlib.pyplot as plt

    sLen = len(scenarioLess)
    row = int(sLen / 2) if sLen % 2 == 0 else int(sLen / 2 + 1)

    for i in range(sLen):
        plt.subplot(row, 2, i+1)
        plt.plot(scenarioGreat[i])
        plt.plot(scenarioLess[i])

    plt.show()

def printIJV(row, dG):
    index = np.where(np.array(dG.i) == row)
    jj = np.array(dG.j)[index]
    vv = np.array(dG.v)[index]
    print(jj)
    print(vv)
