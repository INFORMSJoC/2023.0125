import math
from collections import defaultdict
from itertools import product
from zipfile import ZipFile

import gurobipy as gp
import pandas as pd
from gurobipy import GRB


class AbstractModel:

    def __init__(self):
        self.model = None

    def update(self):
        self.model.update()

    def solve(self, callback=None):
        self.model.optimize(callback=callback)



class DC(AbstractModel):

    def __init__(self, *args, **kwargs):
        # options
        self.options = kwargs['options']
        # sets
        self.N = kwargs['N']
        self.N_G = kwargs['N_G']
        self.N_D = kwargs['N_D']
        self.G = kwargs['G']
        self.D = kwargs['D']
        self.E = kwargs['E']
        self.EA = {(n, m) for (n, m) in kwargs['E']} | {(m, n) for (n, m) in kwargs['E']}
        self.L = kwargs['L']
        self.Omega = kwargs['Omega']
        self.K = kwargs['K']
        self.R = kwargs['R']
        self.KxR = {(k, r) for k in self.R for r in self.R[k]}
        # special sets
        self.NxG = kwargs['NxG']
        self.NxD = kwargs['NxD']
        self.ExL = kwargs['ExL']
        self.G_n = defaultdict(set, kwargs['G_n'])
        self.D_n = defaultdict(set, kwargs['D_n'])
        self.L_nm = defaultdict(set, kwargs['L_nm'])
        self.delta_neg = defaultdict(set, kwargs['delta_neg'])
        self.delta_pos = defaultdict(set, kwargs['delta_pos'])
        self.n_ref = kwargs['n_ref']
        self.k_of_n = kwargs['k_of_n']
        self.NxR = {(n, r) for n in self.N for r in self.R[self.k_of_n[n]]}
        # parameters
        self.p_gen_lo = kwargs['p_gen_lo']
        self.p_gen_hi = kwargs['p_gen_hi']
        self.p_load_hi = kwargs['p_load_hi']
        self.s_flow_hi = kwargs['s_flow_hi']
        self.b = kwargs['b']
        self.probability = kwargs['probability']
        self.xi = defaultdict(lambda: 0, kwargs['xi'])
        self.c = kwargs['c']
        self.f = kwargs['f']
        self.r_hat = kwargs['r_hat']
        self.theta_max = kwargs['options']['theta_max']
        self.theta_delta_max = kwargs['options']['theta_delta_max']
        self.bigM_p_neg = kwargs['bigM']['p_neg']
        self.bigM_p_pos = kwargs['bigM']['p_pos']
        # model
        self.model = model = gp.Model()
        # variables
        _reals_kwargs = {'vtype': GRB.CONTINUOUS, 'lb': -GRB.INFINITY, 'ub': GRB.INFINITY}
        self.p = model.addVars(self.N, self.Omega, **_reals_kwargs, name='p')
        self.p_hat = model.addVars(self.G, self.Omega, **_reals_kwargs, name='p_hat')
        self.p_check = model.addVars(self.G, self.Omega, vtype=GRB.CONTINUOUS, name='p_check')
        self.p_tilde = model.addVars(self.L, self.Omega, **_reals_kwargs, name='p_tilde')
        self.z = model.addVars(self.D, self.Omega, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='z')
        self.zeta = model.addVars(self.NxG, self.Omega, vtype=GRB.BINARY, name='zeta')
        self.theta = model.addVars(self.N, self.Omega, vtype=GRB.CONTINUOUS, lb=-math.pi, ub=math.pi, name='theta')
        self.sin_hat = model.addVars(self.E, self.Omega, **_reals_kwargs, name='sin_hat')
        self.alpha = model.addVars(self.N, self.Omega, vtype=GRB.BINARY, name='alpha')
        self.beta = model.addVars(self.EA, self.Omega, vtype=GRB.BINARY, name='beta')
        self.gamma_under = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma_under')
        self.gamma_over = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma_over')
        self.gamma = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma')
        self.x = model.addVars(self.KxR, vtype=GRB.BINARY, name='x')
        self.chi = model.addVars(self.Omega, vtype=GRB.BINARY, name='chi')
        # constraints
        self.con_resource_hi = model.addConstr(
            sum(self.c[k, r] * self.x[k, r] for k, r in self.KxR) <= self.f,
            name='con_resource_hi'
        )
        self.con_inevitable_damage = model.addConstrs(
            (
                self.x[k, self.r_hat[k]] == 0
                for k in self.K
            ),
            name='con_inevitable_damage'
        )
        self.con_incremental = model.addConstrs(
            (
                self.x[k, r + 1] <= self.x[k, r]
                for (k, r) in self.KxR
                if not r == self.r_hat[k]
            ),
            name='con_incremental'
        )
        self.con_def_alpha_gt = model.addConstrs(
            (
                self.alpha[n, omega] >=\
                sum(1 - self.xi[self.k_of_n[n], r, omega] * (1 - self.x[self.k_of_n[n], r])
                    for r in self.R[self.k_of_n[n]]) -\
                len(self.R[self.k_of_n[n]]) + 1
                for n in self.N
                for omega in self.Omega
            ),
            name='con_def_alpha_gt'
        )
        self.con_def_alpha_lt = model.addConstrs(
            (
                self.alpha[n, omega] <=\
                1 - self.xi[self.k_of_n[n], r, omega] * (1 - self.x[self.k_of_n[n], r])
                for (n, r) in self.NxR
                for omega in self.Omega
            ),
            name='con_def_alpha_lt'
        )
        self.con_def_beta_gt = model.addConstrs(
            (
                self.beta[n, m, omega] >= self.alpha[n, omega] + self.alpha[m, omega] - 1
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_gt'
        )
        self.con_def_beta_lt_f = model.addConstrs(
            (
                self.beta[n, m, omega] <= self.alpha[n, omega]
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_lt_f'
        )
        self.con_def_beta_lt_t = model.addConstrs(
            (
                self.beta[n, m, omega] <= self.alpha[m, omega]
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_lt_t'
        )
        self.con_gen_dispatch = model.addConstrs(
            (
                self.zeta[n, g, omega] == self.alpha[n, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_gen_dispatch'
        )
        self.con_blackout = model.addConstrs(
            (
                self.z[d, omega] <= 1 - self.chi[omega]
                for d in self.D
                for omega in self.Omega
            ),
            name='con_blackout'
        )
        self.con_p_net_flow = model.addConstrs(
            (
                self.p[n, omega] ==\
                sum(self.p_tilde[l, omega]
                    for m in self.delta_neg[n]
                    for l in self.L_nm[m, n]) -\
                sum(self.p_tilde[l, omega]
                    for m in self.delta_pos[n]
                    for l in self.L_nm[n, m])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_p_net_flow'
        )
        self.con_p_net_injection = model.addConstrs(
            (
                self.p[n, omega] ==\
                sum(self.p_hat[g, omega] - self.p_check[g, omega]
                    for g in self.G_n[n]) -\
                sum(self.p_load_hi[d] * self.z[d, omega]
                    for d in self.D_n[n])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_p_net_injection'
        )
        self.con_p_hat_p_check = model.addConstrs(
            (
                self.p_check[g, omega] <= self.p_hat[g, omega]
                for g in self.G
                for omega in self.Omega
            ),
            name='con_p_hat_p_check'
        )
        self.con_p_gen_lo = model.addConstrs(
            (
                self.p_gen_lo[g] * (self.zeta[n, g, omega] - self.chi[omega]) <= self.p_hat[g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_p_gen_lo'
        )
        self.con_p_gen_hi = model.addConstrs(
            (
                self.p_hat[g, omega] <= self.p_gen_hi[g] * self.zeta[n, g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_p_gen_hi'
        )
        self.con_p_flow_lo = model.addConstrs(
            (
                -self.s_flow_hi[l] * self.beta[n, m, omega] <= self.p_tilde[l, omega]
                for (n, m, l) in self.ExL
                for omega in self.Omega
            ),
            name='con_p_flow_lo'
        )
        self.con_p_flow_hi = model.addConstrs(
            (
                self.p_tilde[l, omega] <= self.s_flow_hi[l] * self.beta[n, m, omega]
                for (n, m, l) in self.ExL
                for omega in self.Omega
            ),
            name='con_p_flow_hi'
        )
        self.con_ohms_law_gt = model.addConstrs(
            (
                -self.b[l] * self.sin_hat[n, m, omega] >=\
                self.p_tilde[l, omega] + self.bigM_p_neg[n, m, l] * (1 - self.beta[n, m, omega])
                for (n, m, l) in self.ExL
                for omega in self.Omega
            ),
            name='con_ohms_law_gt'
        )
        self.con_ohms_law_lt = model.addConstrs(
            (
                -self.b[l] * self.sin_hat[n, m, omega] <=\
                self.p_tilde[l, omega] + self.bigM_p_pos[n, m, l] * (1 - self.beta[n, m, omega])
                for (n, m, l) in self.ExL
                for omega in self.Omega
            ),
            name='con_ohms_law_lt'
        )
        # linear approximation of sine
        self.con_sin_hat = model.addConstrs(
            (
                self.sin_hat[n, m, omega] == self.theta[n, omega] - self.theta[m, omega]
                for (n, m) in self.E
                for omega in self.Omega
            ),
            name='con_sin_hat'
        )
        self.con_sin_hat_lo = model.addConstrs(
            (
                -2 * (1 - self.beta[n, m, omega]) * self.theta_max\
                - self.beta[n, m, omega] * self.theta_delta_max <=\
                self.sin_hat[n, m, omega]
                for (n, m) in self.E
                for omega in self.Omega
            ),
            name='con_sin_hat_lo'
        )
        self.con_sin_hat_hi = model.addConstrs(
            (
                self.sin_hat[n, m, omega] <=\
                2 * (1 - self.beta[n, m, omega]) * self.theta_max\
                + self.beta[n, m, omega] * self.theta_delta_max
                for (n, m) in self.E
                for omega in self.Omega
            ),
            name='con_sin_hat_hi'
        )
        self.con_ref_phase_angle = model.addConstrs(
            (
                self.theta[self.n_ref, omega] == 0
                for omega in self.Omega
            ),
            name='con_ref_phase_angle'
        )
        self.con_def_gamma_under = model.addConstrs(
            (
                self.gamma_under[omega] ==\
                sum(self.p_load_hi[d] * (1 - self.z[d, omega]) for d in self.D)
                for omega in self.Omega
            ),
            name='con_def_gamma_under'
        )
        self.con_def_gamma_over = model.addConstrs(
            (
                self.gamma_over[omega] ==\
                sum(self.p_check[g, omega] for g in self.G)
                for omega in self.Omega
            ),
            name='con_def_gamma_over'
        )
        self.con_def_gamma = model.addConstrs(
            (
                self.gamma[omega] == self.gamma_under[omega] + self.gamma_over[omega]
                for omega in self.Omega
            ),
            name='con_def_gamma'
        )
        if self.options['approach'] == 'stochastic':
            self.obj_min_expected_shed =\
                sum(self.probability[omega] * self.gamma[omega]
                    for omega in self.Omega)
            model.setObjective(self.obj_min_expected_shed, sense=GRB.MINIMIZE)
        elif self.options['approach'] == 'robust':
            self.gamma_max = model.addVar(vtype=GRB.CONTINUOUS, name='gamma_max')
            self.con_gamma_max = model.addConstrs(
                self.gamma_max >= self.gamma[omega]
                for omega in self.Omega
            )
            self.obj_min_max_shed = self.gamma_max
            model.setObjective(self.obj_min_max_shed, sense=GRB.MINIMIZE)
        else:
            raise ValueError("Option 'approach' must be one of 'stochastic' or 'robust'.")


class LPAC(AbstractModel):

    def __init__(self, *args, **kwargs):
        # options
        self.options = kwargs['options']
        # sets
        self.N = kwargs['N']
        self.N_G = kwargs['N_G']
        self.N_D = kwargs['N_D']
        self.G = kwargs['G']
        self.D = kwargs['D']
        self.E = kwargs['E']
        self.EA = {(n, m) for (n, m) in kwargs['E']} | {(m, n) for (n, m) in kwargs['E']}
        self.L = kwargs['L']
        self.O = {'f', 'b'}
        self.Omega = kwargs['Omega']
        self.K = kwargs['K']
        self.R = kwargs['R']
        self.KxR = {(k, r) for k in self.R for r in self.R[k]}
        # special sets
        self.NxG = kwargs['NxG']
        self.NxD = kwargs['NxD']
        self.EAxLxO = {(n, m, l, 'f') for (n, m, l) in kwargs['ExL']}\
                    | {(m, n, l, 'b') for (n, m, l) in kwargs['ExL']}
        self.G_n = defaultdict(set, kwargs['G_n'])
        self.D_n = defaultdict(set, kwargs['D_n'])
        self.L_nm = defaultdict(set, kwargs['L_nm'])
        self.LxO_nm = defaultdict(set)
        for (n, m) in kwargs['E']:
            for l in kwargs['L_nm'][n, m]:
                self.LxO_nm[n, m].add((l, 'f'))
                self.LxO_nm[m, n].add((l, 'b'))
        self.delta = defaultdict(set)
        for n in kwargs['delta_neg']:
            self.delta[n] |= kwargs['delta_neg'][n]
        for n in kwargs['delta_pos']:
            self.delta[n] |= kwargs['delta_pos'][n]
        self.n_ref = kwargs['n_ref']
        self.k_of_n = kwargs['k_of_n']
        self.NxR = {(n, r) for n in self.N for r in self.R[self.k_of_n[n]]}
        # parameters
        self.p_gen_lo = kwargs['p_gen_lo']
        self.p_gen_hi = kwargs['p_gen_hi']
        self.q_gen_lo = kwargs['q_gen_lo']
        self.q_gen_hi = kwargs['q_gen_hi']
        self.p_load_hi = kwargs['p_load_hi']
        self.q_load_hi = kwargs['q_load_hi']
        self.s_flow_hi = kwargs['s_flow_hi']
        self.b = kwargs['b']
        self.g = kwargs['g']
        self.v = kwargs['v']
        self.v_lo = kwargs['v_lo']
        self.v_hi = kwargs['v_hi']
        self.probability = kwargs['probability']
        self.xi = defaultdict(lambda: 0, kwargs['xi'])
        self.c = kwargs['c']
        self.f = kwargs['f']
        self.r_hat = kwargs['r_hat']
        self.theta_max = kwargs['options']['theta_max']
        self.theta_delta_max = kwargs['options']['theta_delta_max']
        self.bigM_p_neg = kwargs['bigM']['p_neg']
        self.bigM_p_pos = kwargs['bigM']['p_pos']
        self.bigM_q_neg = kwargs['bigM']['q_neg']
        self.bigM_q_pos = kwargs['bigM']['q_pos']
        # model
        self.model = model = gp.Model()
        # variables
        _reals_kwargs = {'vtype': GRB.CONTINUOUS, 'lb': -GRB.INFINITY, 'ub': GRB.INFINITY}
        self.p = model.addVars(self.N, self.Omega, **_reals_kwargs, name='p')
        self.q = model.addVars(self.N, self.Omega, **_reals_kwargs, name='q')
        self.p_hat = model.addVars(self.G, self.Omega, **_reals_kwargs, name='p_hat')
        self.q_hat = model.addVars(self.G, self.Omega, **_reals_kwargs, name='q_hat')
        self.p_check = model.addVars(self.G, self.Omega, vtype=GRB.CONTINUOUS, name='p_check')
        self.p_tilde = model.addVars(self.L, self.O, self.Omega, **_reals_kwargs, name='p_tilde')
        self.q_tilde = model.addVars(self.L, self.O, self.Omega, **_reals_kwargs, name='q_tilde')
        self.z = model.addVars(self.D, self.Omega, vtype=GRB.CONTINUOUS, lb=0.0, ub=1.0, name='z')
        self.zeta = model.addVars(self.NxG, self.Omega, vtype=GRB.BINARY, name='zeta')
        self.theta = model.addVars(self.N, self.Omega, vtype=GRB.CONTINUOUS, lb=-math.pi, ub=math.pi, name='theta')
        self.phi = model.addVars(self.N, self.Omega, **_reals_kwargs, name='phi')
        self.sin_hat = model.addVars(self.EA, self.Omega, **_reals_kwargs, name='sin_hat')
        self.cos_hat = model.addVars(self.EA, self.Omega, vtype=GRB.CONTINUOUS, lb=math.cos(self.theta_delta_max), ub=1, name='cos_hat')
        self.alpha = model.addVars(self.N, self.Omega, vtype=GRB.BINARY, name='alpha')
        self.beta = model.addVars(self.EA, self.Omega, vtype=GRB.BINARY, name='beta')
        self.gamma_under = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma_under')
        self.gamma_over = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma_over')
        self.gamma = model.addVars(self.Omega, vtype=GRB.CONTINUOUS, name='gamma')
        self.x = model.addVars(self.KxR, vtype=GRB.BINARY, name='x')
        self.chi = model.addVars(self.Omega, vtype=GRB.BINARY, name='chi')
        # constraints
        self.con_resource_hi = model.addConstr(
            sum(self.c[k, r] * self.x[k, r] for k, r in self.KxR) <= self.f,
            name='con_resource_hi'
        )
        self.con_inevitable_damage = model.addConstrs(
            (
                self.x[k, self.r_hat[k]] == 0
                for k in self.K
            ),
            name='con_inevitable_damage'
        )
        self.con_incremental = model.addConstrs(
            (
                self.x[k, r + 1] <= self.x[k, r]
                for (k, r) in self.KxR
                if not r == self.r_hat[k]
            ),
            name='con_incremental'
        )
        self.con_def_alpha_gt = model.addConstrs(
            (
                self.alpha[n, omega] >=\
                sum(1 - self.xi[self.k_of_n[n], r, omega] * (1 - self.x[self.k_of_n[n], r])
                    for r in self.R[self.k_of_n[n]]) -\
                len(self.R[self.k_of_n[n]]) + 1
                for n in self.N
                for omega in self.Omega
            ),
            name='con_def_alpha_gt'
        )
        self.con_def_alpha_lt = model.addConstrs(
            (
                self.alpha[n, omega] <=\
                1 - self.xi[self.k_of_n[n], r, omega] * (1 - self.x[self.k_of_n[n], r])
                for (n, r) in self.NxR
                for omega in self.Omega
            ),
            name='con_def_alpha_lt'
        )
        self.con_def_beta_gt = model.addConstrs(
            (
                self.beta[n, m, omega] >= self.alpha[n, omega] + self.alpha[m, omega] - 1
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_gt'
        )
        self.con_def_beta_lt_f = model.addConstrs(
            (
                self.beta[n, m, omega] <= self.alpha[n, omega]
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_lt_f'
        )
        self.con_def_beta_lt_t = model.addConstrs(
            (
                self.beta[n, m, omega] <= self.alpha[m, omega]
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_def_beta_lt_t'
        )
        self.con_gen_dispatch = model.addConstrs(
            (
                self.zeta[n, g, omega] == self.alpha[n, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_gen_dispatch'
        )
        self.con_blackout = model.addConstrs(
            (
                self.z[d, omega] <= 1 - self.chi[omega]
                for d in self.D
                for omega in self.Omega
            ),
            name='con_blackout'
        )
        # active power constraints
        self.con_p_net_flow = model.addConstrs(
            (
                self.p[n, omega] ==\
                sum(self.p_tilde[l, o, omega]
                    for m in self.delta[n]
                    for (l, o) in self.LxO_nm[n, m])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_p_net_flow'
        )
        self.con_p_net_injection = model.addConstrs(
            (
                self.p[n, omega] ==\
                sum(self.p_hat[g, omega] - self.p_check[g, omega]
                    for g in self.G_n[n]) -\
                sum(self.p_load_hi[d] * self.z[d, omega]
                    for d in self.D_n[n])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_p_net_injection'
        )
        self.con_p_hat_p_check = model.addConstrs(
            (
                self.p_check[g, omega] <= self.p_hat[g, omega]
                for g in self.G
                for omega in self.Omega
            ),
            name='con_p_hat_p_check'
        )
        self.con_p_gen_lo = model.addConstrs(
            (
                self.p_gen_lo[g] * (self.zeta[n, g, omega] - self.chi[omega]) <= self.p_hat[g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_p_gen_lo'
        )
        self.con_p_gen_hi = model.addConstrs(
            (
                self.p_hat[g, omega] <= self.p_gen_hi[g] * self.zeta[n, g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_p_gen_hi'
        )
        self.con_p_ohms_law_gt = model.addConstrs(
            (
                self.p_tilde[l, o, omega] >=\
                    self.v[n] * self.g[l] * (self.v[m] - self.v[n]) * self.chi[omega] +\
                    self.v[n] ** 2 * self.g[l] -\
                    self.v[n] * self.v[m] * self.g[l] * self.cos_hat[n, m, omega] +\
                    self.v[n] * self.v[m] * self.b[l] * self.sin_hat[n, m, omega] +\
                    self.bigM_p_neg[n, m, l, o] * (1 - self.beta[n, m, omega])
                for n, m, l, o in self.EAxLxO
                for omega in self.Omega
            ),
            name='con_p_ohms_law_gt'
        )
        self.con_p_ohms_law_lt = model.addConstrs(
            (
                self.p_tilde[l, o, omega] <=\
                    self.v[n] * self.g[l] * (self.v[m] - self.v[n]) * self.chi[omega] +\
                    self.v[n] ** 2 * self.g[l] -\
                    self.v[n] * self.v[m] * self.g[l] * self.cos_hat[n, m, omega] +\
                    self.v[n] * self.v[m] * self.b[l] * self.sin_hat[n, m, omega] +\
                    self.bigM_p_pos[n, m, l, o] * (1 - self.beta[n, m, omega])
                for n, m, l, o in self.EAxLxO
                for omega in self.Omega
            ),
            name='con_p_ohms_law_lt'
        )
        # reactive power constraints
        self.con_q_net_flow = model.addConstrs(
            (
                self.q[n, omega] ==\
                sum(self.q_tilde[l, o, omega]
                    for m in self.delta[n]
                    for (l, o) in self.LxO_nm[n, m])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_q_net_flow'
        )
        self.con_q_net_injection = model.addConstrs(
            (
                self.q[n, omega] ==\
                sum(self.q_hat[g, omega]
                    for g in self.G_n[n]) -\
                sum(self.q_load_hi[d] * self.z[d, omega]
                    for d in self.D_n[n])
                for n in self.N
                for omega in self.Omega
            ),
            name='con_q_net_injection'
        )
        self.con_q_gen_lo = model.addConstrs(
            (
                self.q_gen_lo[g] * self.zeta[n, g, omega] <= self.q_hat[g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_q_gen_lo'
        )
        self.con_q_gen_hi = model.addConstrs(
            (
                self.q_hat[g, omega] <= self.q_gen_hi[g] * self.zeta[n, g, omega]
                for (n, g) in self.NxG
                for omega in self.Omega
            ),
            name='con_q_gen_hi'
        )
        self.con_q_ohms_law_gt = model.addConstrs(
            (
                self.q_tilde[l, o, omega] >=\
                    self.v[n] * self.b[l] * (self.v[n] - self.v[m]) * self.chi[omega] -\
                    self.v[n] ** 2 * self.b[l] -\
                    self.v[n] * self.v[m] * self.g[l] * self.sin_hat[n, m, omega] +\
                    self.v[n] * self.v[m] * self.b[l] * self.cos_hat[n, m, omega] -\
                    self.v[n] * self.b[l] * (self.phi[n, omega] - self.phi[m, omega]) -\
                    (self.v[n] - self.v[m]) * self.b[l] * self.phi[n, omega] +\
                    self.bigM_q_neg[n, m, l, o] * (1 - self.beta[n, m, omega])
                for n, m, l, o in self.EAxLxO
                for omega in self.Omega
            ),
            name='con_q_ohms_law_gt'
        )
        self.con_q_ohms_law_lt = model.addConstrs(
            (
                self.q_tilde[l, o, omega] <=\
                    self.v[n] * self.b[l] * (self.v[n] - self.v[m]) * self.chi[omega] -\
                    self.v[n] ** 2 * self.b[l] -\
                    self.v[n] * self.v[m] * self.g[l] * self.sin_hat[n, m, omega] +\
                    self.v[n] * self.v[m] * self.b[l] * self.cos_hat[n, m, omega] -\
                    self.v[n] * self.b[l] * (self.phi[n, omega] - self.phi[m, omega]) -\
                    (self.v[n] - self.v[m]) * self.b[l] * self.phi[n, omega] +\
                    self.bigM_q_pos[n, m, l, o] * (1 - self.beta[n, m, omega])
                for n, m, l, o in self.EAxLxO
                for omega in self.Omega
            ),
            name='con_q_ohms_law_lt'
        )
        # linear approximation of sine
        self.con_sin_hat = model.addConstrs(
            (
                self.sin_hat[n, m, omega] == self.theta[n, omega] - self.theta[m, omega]
                for (n, m) in self.EA
                for omega in self.Omega
            ),
            name='con_sin_hat'
        )
        self.con_sin_hat_lo = model.addConstrs(
            (
                -2 * (1 - self.beta[n, m, omega]) * self.theta_max\
                - self.beta[n, m, omega] * self.theta_delta_max <=\
                self.sin_hat[n, m, omega]
                for (n, m) in self.E
                for omega in self.Omega
            ),
            name='con_sin_hat_lo'
        )
        self.con_sin_hat_hi = model.addConstrs(
            (
                self.sin_hat[n, m, omega] <=\
                2 * (1 - self.beta[n, m, omega]) * self.theta_max\
                + self.beta[n, m, omega] * self.theta_delta_max
                for (n, m) in self.E
                for omega in self.Omega
            ),
            name='con_sin_hat_hi'
        )
        # cosine polyhedral relaxation
        if self.options['cosine_model'] == 'constant':
            self.con_cos_hat = model.addConstrs(
                (
                    self.cos_hat[n, m, omega] == 1
                    for (n, m) in self.EA
                    for omega in self.Omega
                ),
                name='con_cos_hat'
            )
        elif self.options['cosine_model'] == 'linear':
            f = lambda theta, theta_hat: (theta_hat - theta) * math.sin(theta_hat) + math.cos(theta_hat)
            self.con_cos_hat = model.addConstrs(
                (
                    self.cos_hat[n, m, omega] <=\
                    f(self.theta[n, omega] - self.theta[m, omega], t)\
                    + (1 - self.beta[n, m, omega]) * (1 - min(f(-2 * self.theta_max, t), f(2 * self.theta_max, t)))
                    for (n, m) in self.EA
                    for t in self.options['Theta_cos']
                    for omega in self.Omega
                ),
                name='con_cos_hat'
            )
        elif self.options['cosine_model'] == 'quadratic':
            g = lambda theta: 1 - (1 - math.cos(self.theta_delta_max)) * (theta / self.theta_delta_max) ** 2
            self.con_cos_hat = model.addConstrs(
                (
                    self.cos_hat[n, m, omega] <=\
                    g(self.theta[n, omega] - self.theta[m, omega])\
                    + (1 - self.beta[n, m, omega]) * (1 - min(g(-2 * self.theta_max), g(2 * self.theta_max)))
                    for (n, m) in self.EA
                    for omega in self.Omega
                ),
                name='con_cos_hat'
            )
        else:
            raise ValueError("Option 'cosine_model' must be one of 'constant', 'linear', or 'quadratic'.")
        # thermal limit relaxation
        if self.options['flow_limit_model'] == 'outer':
            self.con_flow_limit = model.addConstrs(
                (
                    math.cos(t) * self.p_tilde[l, o, omega] +\
                    math.sin(t) * self.q_tilde[l, o, omega] <=\
                        self.s_flow_hi[l] * self.beta[n, m, omega]
                    for t in self.options['Theta_flow_limit']
                    for (n, m, l, o) in self.EAxLxO
                    for omega in self.Omega
                ),
                name='con_flow_limit'
            )
        elif self.options['flow_limit_model'] == 'inner':
            # this may not be correct
            self.con_flow_limit = model.addConstrs(
                (
                    math.cos(t + math.pi / len(self.options['Theta_flow_limit'])) / math.cos(math.pi / len(self.options['Theta_flow_limit'])) * self.p_tilde[l, o, omega] -\
                    math.sin(t + math.pi / len(self.options['Theta_flow_limit'])) / math.cos(math.pi / len(self.options['Theta_flow_limit'])) * self.p_tilde[l, o, omega] <=\
                        self.s_flow_hi[l] * self.beta[n, m, omega]
                    for t in self.options['Theta_flow_limit']
                    for (n, m, l, o) in self.EAxLxO
                    for omega in self.Omega
                ),
                name='con_flow_limit'
            )
        elif self.options['flow_limit_model'] == 'exact':
            self.con_flow_limit = model.addConstrs(
                (
                    self.p_tilde[l, o, omega] ** 2 + self.q_tilde[l, o, omega] ** 2 <=\
                        self.s_flow_hi[l] ** 2 * self.beta[n, m, omega]
                    for (n, m, l, o) in self.EAxLxO
                    for omega in self.Omega
                ),
                name='con_flow_limit'
            )
        else:
            raise ValueError("Option 'flow_limit_model' must be one of 'outer', 'inner', or 'exact'.")
        # other
        self.con_voltage_lo = model.addConstrs(
            (
                self.v_lo[n] <= self.v[n] + self.phi[n, omega]
                for n in self.N
                for omega in self.Omega
            ),
            name='con_voltage_lo'
        )
        self.con_voltage_hi = model.addConstrs(
            (
                self.v[n] + self.phi[n, omega] <= self.v_hi[n]
                for n in self.N
                for omega in self.Omega
            ),
            name='con_voltage_hi'
        )
        self.con_ref_phase_angle = model.addConstrs(
            (
                self.theta[self.n_ref, omega] == 0
                for omega in self.Omega
            ),
            name='con_ref_phase_angle'
        )
        self.con_ref_voltage = model.addConstrs(
            (
                self.phi[self.n_ref, omega] == 0
                for omega in self.Omega
            ),
            name='con_ref_voltage'
        )
        self.con_def_gamma_under = model.addConstrs(
            (
                self.gamma_under[omega] ==\
                sum(self.p_load_hi[d] * (1 - self.z[d, omega]) for d in self.D)
                for omega in self.Omega
            ),
            name='con_def_gamma_under'
        )
        self.con_def_gamma_over = model.addConstrs(
            (
                self.gamma_over[omega] ==\
                sum(self.p_check[g, omega] for g in self.G)
                for omega in self.Omega
            ),
            name='con_def_gamma_over'
        )
        self.con_def_gamma = model.addConstrs(
            (
                self.gamma[omega] == self.gamma_under[omega] + self.gamma_over[omega]
                for omega in self.Omega
            ),
            name='con_def_gamma'
        )
        if self.options['approach'] == 'stochastic':
            self.obj_min_expected_shed =\
                sum(self.probability[omega] * self.gamma[omega]
                    for omega in self.Omega)
            model.setObjective(self.obj_min_expected_shed, sense=GRB.MINIMIZE)
        elif self.options['approach'] == 'robust':
            self.gamma_max = model.addVar(vtype=GRB.CONTINUOUS, name='gamma_max')
            self.con_gamma_max = model.addConstrs(
                self.gamma_max >= self.gamma[omega]
                for omega in self.Omega
            )
            self.obj_min_max_shed = self.gamma_max
            model.setObjective(self.obj_min_max_shed, sense=GRB.MINIMIZE)
        else:
            raise ValueError("Option 'approach' must be one of 'stochastic' or 'robust'.")


class Solution:

    def __init__(self, solution):
        self._solution = solution

    def __getitem__(self, variable):
        return self._solution[variable]

    def keys(self):
        return self._solution.keys()

    @classmethod
    def from_solved_instance(cls, instance):
        solution = {}
        for key1, val1 in vars(instance).items():
            if isinstance(val1, gp.tupledict):
                key2 = next(iter(val1))
                if isinstance(val1[key2], gp.Var):
                    solution[key1] = pd.Series({key2: val2.X for key2, val2 in val1.items()})
        solution['ObjVal'] = instance.model.ObjVal
        solution['ObjBound'] = instance.model.ObjBound
        return cls(solution)

    @classmethod
    def from_zip(cls, filename):
        solution = {}
        with ZipFile(filename) as zfh:
            for obj in zfh.filelist:
                if obj.filename.endswith('.csv'):
                    df = pd.read_csv(zfh.open(obj.filename), header=None)
                    df.set_index(list(df.columns[:-1]), inplace=True)
                    df.index.names = [None] * len(df.index.names)
                    ser = pd.Series(df[df.columns[0]])
                    ser.name = None
                    solution[obj.filename.split('.')[0]] = ser
                elif obj.filename.endswith('.txt'):
                    with zfh.open(obj.filename) as fh:
                        solution[obj.filename.split('.')[0]] = float(fh.read().strip())
                else:
                    continue
        return cls(solution)

    def to_zip(self, filename, variables=None):
        if variables is None:
            variables = list(self._solution)
        with ZipFile(filename, 'w') as zfh:
            for variable in variables:
                if variable in self._solution:
                    if not type(self._solution[variable]) is pd.Series:
                        zfh.writestr(f'{variable}.txt', str(self._solution[variable]))
                    else:
                        zfh.writestr(f'{variable}.csv', self._solution[variable].to_csv(header=False))
