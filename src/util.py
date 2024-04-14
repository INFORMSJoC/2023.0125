import os
from collections import defaultdict
from copy import deepcopy
from itertools import combinations
from zipfile import ZipFile

import numpy as np
import pandas as pd
import pyomo.environ as pe
import yaml

from config import EMAIL, SPECS_FLOOD_YAMLS, CONSOLIDATED_RESULTS_DIR


def in_notebook():
    try:
        from IPython import get_ipython
        if 'IPKernelApp' not in get_ipython().config:  # pragma: no cover
            return False
    except ImportError:
        return False
    except AttributeError:
        return False
    return True


def max_budget_by_omega(c, xi, r_hat, **specs):
    maxlevels = {(k, omega): r for (k, r, omega) in sorted(xi.keys())}
    Omega = {omega for (_, _, omega) in xi.keys()}
    newlevels = {(k, omega, rp)
                 for (k, omega), r in maxlevels.items()
                 for rp in range(1, r + 1)
                 if r < r_hat[k]}
    c_star = {omega: sum(c[k, r] for (k, omegap, r) in newlevels if omegap == omega)
              for omega in Omega}
    return c_star


class Grid:

    def __init__(self, specs):
        self.xi = dict()
        for omega in specs['Omega']:
            self.xi[omega] = dict()
            for k in specs['K']:
                for r in specs['R'][k]:
                    if (k, r, omega) in specs['xi']:
                        self.xi[omega][k, r] = 1
                    else:
                        self.xi[omega][k, r] = 0
        self.N_k = defaultdict(set)
        for n, k in specs['k_of_n'].items():
            self.N_k[k].add(n)
        self.delta_pos = defaultdict(set, deepcopy(specs['delta_pos']))
        self.delta_neg = defaultdict(set, deepcopy(specs['delta_neg']))
        for key, val in specs.items():
            if key not in ['xi', 'delta_pos', 'delta_neg']:
                setattr(self, key, val)

    def get_cost(self, decision):
        return sum(self.c[action] for action in decision)


class ResilienceState:

    def __init__(self, grid):
        self.grid = grid
        self.x = dict()
        self.alpha = {omega: dict() for omega in self.grid.Omega}
        self.beta = {omega: dict() for omega in self.grid.Omega}
        self._outed_load = None
        self._outed_flow = None
        self.reset()

    def _expand_combo_into_decision(self, combo):
        return tuple((k, rp)
                     for (k, r) in combo
                     for rp in sorted(self.grid.R[k], reverse=True)
                     if rp <= r and self.x[k][rp] == 0)

    def reset(self):
        self._outed_load = {omega: 0 for omega in self.grid.Omega}
        self._outed_flow = {omega: 0 for omega in self.grid.Omega}
        for k in self.grid.K:
            for r in self.grid.R[k]:
                if k not in self.x:
                    self.x[k] = dict()
                self.x[k][r] = 0
        for omega in self.grid.Omega:
            for n in self.grid.N:
                k = self.grid.k_of_n[n]
                if k in self.x:
                    update = 1
                    for r in self.grid.R[k]:
                        update *= (1 - self.grid.xi[omega][k, r] * (1 - self.x[k][r]))
                    self.alpha[omega][n] = update
                else:
                    self.alpha[omega][n] = 1
                if self.alpha[omega][n] == 0:
                    self._outed_load[omega] += sum(self.grid.p_load_hi[d]
                                                   for d in self.grid.D_n.get(n, set()))
        for omega in self.grid.Omega:
            for n, m in self.grid.E:
                self.beta[omega][n, m] = self.alpha[omega][n] * self.alpha[omega][m]
                if self.beta[omega][n, m] == 0:
                    self._outed_flow[omega] += sum(self.grid.s_flow_hi[l]
                                                   for l in self.grid.L_nm.get((n, m), set())
                                                   if self.grid.s_flow_hi[l] < 10000)

    def get_cost(self, decision):
        return self.grid.get_cost(decision)

    def get_stochastic_benefit(self, decision, weight_load=1, weight_flow=0):
        x_diff = self.get_x_diff(decision)
        alpha_diff = self.get_alpha_diff(x_diff)
        if weight_load != 0:
            benefit_load = sum(self.grid.probability[omega] * self.grid.p_load_hi[d]
                               for omega in alpha_diff
                               for n in alpha_diff[omega]
                               for d in self.grid.D_n.get(n, set()))
        else:
            benefit_load = 0
        if weight_flow != 0:
            beta_diff = self.get_beta_diff(alpha_diff)
            benefit_flow = sum(self.grid.probability[omega] * self.grid.s_flow_hi[l]
                               for omega in beta_diff
                               for n, m in beta_diff[omega]
                               for l in self.grid.L_nm.get((n, m), set())
                               if self.grid.s_flow_hi[l] < 10000)
        else:
            benefit_flow = 0
        return weight_load * benefit_load + weight_flow * benefit_flow

    def get_robust_benefit(self, decision, weight_load=1, weight_flow=0):
        def omega_star_key(omega):
            return weight_load * self._outed_load[omega] + weight_flow * self._outed_flow[omega]
        omega_star = max(self.grid.Omega, key=omega_star_key)
        x_diff = self.get_x_diff(decision)
        alpha_diff = self.get_alpha_diff(x_diff)
        if weight_load != 0:
            benefit_load = sum(self.grid.p_load_hi[d]
                               for n in alpha_diff[omega_star]
                               for d in self.grid.D_n.get(n, set()))
        else:
            benefit_load = 0
        if weight_flow != 0:
            beta_diff = self.get_beta_diff(alpha_diff)
            benefit_flow = sum(self.grid.s_flow_hi[l]
                               for n, m in beta_diff[omega_star]
                               for l in self.grid.L_nm.get((n, m), set())
                               if self.grid.s_flow_hi[l] < 10000)
        else:
            benefit_flow = 0
        return weight_load * benefit_load + weight_flow * benefit_flow

    def get_stochastic_benefit_to_cost_ratio(self, decision, weight_load=1, weight_flow=0):
        benefit = self.get_stochastic_benefit(decision, weight_load=weight_load, weight_flow=weight_flow)
        cost = self.get_cost(decision)
        return benefit / cost

    def get_robust_benefit_to_cost_ratio(self, decision, weight_load=1, weight_flow=0):
        benefit = self.get_stochastic_benefit(decision, weight_load=weight_load, weight_flow=weight_flow)
        cost = self.get_cost(decision)
        return benefit / cost

    def get_x_diff(self, decision):
        x_diff = defaultdict(dict)
        for (k, r) in decision:
            if self.x[k][r] == 1:
                raise ValueError(f'Action ({k}, {r}) already implemented.')
            x_diff[k][r] = 1
        return x_diff

    def get_alpha_diff(self, x_diff):
        alpha_diff = defaultdict(dict)
        for omega in self.grid.Omega:
            for k in x_diff:
                for n in self.grid.N_k[k]:
                    update = 1
                    for r in self.grid.R[k]:
                        update *= (1 - self.grid.xi[omega][k, r] * (1 - {**self.x[k], **x_diff[k]}[r]))
                    if self.alpha[omega][n] != update:
                        alpha_diff[omega][n] = update
        return alpha_diff

    def get_beta_diff(self, alpha_diff):
        beta_diff = defaultdict(dict)
        for omega in alpha_diff:
            for n in alpha_diff[omega]:
                for m in self.grid.delta_pos[n] | self.grid.delta_neg[n]:
                    update = 1
                    update *= {**self.alpha[omega], **alpha_diff[omega]}[n]
                    update *= {**self.alpha[omega], **alpha_diff[omega]}[m]
                    no, mo = (n, m) if n < m else (m, n)
                    if self.beta[omega][no, mo] != update:
                        beta_diff[omega][no, mo] = update
        return beta_diff

    def adopt_decision(self, decision):
        x_diff = self.get_x_diff(decision)
        for k in x_diff:
            self.x[k].update(x_diff[k])
        alpha_diff = self.get_alpha_diff(x_diff)
        for omega in alpha_diff:
            self.alpha[omega].update(alpha_diff[omega])
            for n in alpha_diff[omega]:
                self._outed_load[omega] -= sum(self.grid.p_load_hi[d]
                                               for d in self.grid.D_n.get(n, set()))
        beta_diff = self.get_beta_diff(alpha_diff)
        for omega in beta_diff:
            self.beta[omega].update(beta_diff[omega])
            for (n, m) in beta_diff[omega]:
                self._outed_flow[omega] -= sum(self.grid.s_flow_hi[l]
                                               for l in self.grid.L_nm.get((n, m), set())
                                               if self.grid.s_flow_hi[l] < 10000)

    def get_remaining_actions(self):
        for omega in self.grid.xi:
            for (k, r) in self.grid.xi[omega]:
                condition1 = self.grid.xi[omega][k, r] == 1
                condition2 = self.grid.xi[omega][k, self.grid.r_hat[k]] == 0
                condition3 = self.x[k][r] == 0
                if condition1 and condition2 and condition3:
                    yield (k, r)

    def get_cost_feasible_decisions(self, remaining_budget, degree=1):
        for d in range(1, degree + 1):
            for combo in combinations(self.get_remaining_actions(), r=d):
                decision = self._expand_combo_into_decision(combo)
                if self.get_cost(decision) <= remaining_budget:
                    yield decision

    def get_stochastic_greedy_decision(self, budget, degree=1, weight_load=1, weight_flow=0, overspend=False):
        # It may occur that none of the remaining cost-feasible actions provides
        # a benefit. For example, consider what happens if `weight_flow=0` and
        # no still-affected substations have a load. As a pathological example,
        # consider also the case in which the only remaining decision is to
        # enable a single-bus substation with no load and whose neighbors are
        # all down due to inexorable flooding in every scenario such that no
        # branch capacity may be enabled a result of protecting the substation.
        # Setting `overspend = True` will allow the action of protecting a
        # zero-benefit substation. However, note that hardening to a
        # substation's inexorable flood level is not made possible by setting
        # `overspend = True`.
        decision_star = None
        ratio_star = -1 if overspend else 0
        for decision in self.get_cost_feasible_decisions(budget, degree=degree):
            ratio = self.get_stochastic_benefit_to_cost_ratio(decision, weight_load=weight_load, weight_flow=weight_flow)
            if ratio > ratio_star:
                decision_star = decision
                ratio_star = ratio
        return decision_star

    def get_robust_greedy_decision(self, budget, degree=1, weight_load=1, weight_flow=0, overspend=False):
        decision_star = None
        ratio_star = -1 if overspend else 0
        for decision in self.get_cost_feasible_decisions(budget, degree=degree):
            ratio = self.get_robust_benefit_to_cost_ratio(decision, weight_load=weight_load, weight_flow=weight_flow)
            if ratio > ratio_star:
                decision_star = decision
                ratio_star = ratio
        return decision_star

    def adopt_stochastic_greedy_solution(self, budget, degree=1, weight_load=1, weight_flow=0, overspend=False):
        self.reset()
        remaining_budget = budget
        solution = list()
        while remaining_budget > 0:
            decision = self.get_stochastic_greedy_decision(remaining_budget,
                                                           degree=degree,
                                                           weight_load=weight_load,
                                                           weight_flow=weight_flow,
                                                           overspend=overspend)
            if decision is None:
                break
            cost = self.get_cost(decision)
            solution.append(decision)
            self.adopt_decision(decision)
            remaining_budget -= cost
        return solution

    def adopt_robust_greedy_solution(self, budget, degree=1, weight_load=1, weight_flow=0, overspend=False):
        self.reset()
        remaining_budget = budget
        solution = list()
        while remaining_budget > 0:
            decision = self.get_robust_greedy_decision(remaining_budget,
                                                       degree=degree,
                                                       weight_load=weight_load,
                                                       weight_flow=weight_flow,
                                                       overspend=overspend)
            if decision is None:
                break
            cost = self.get_cost(decision)
            solution.append(decision)
            self.adopt_decision(decision)
            remaining_budget -= cost
        return solution


def get_Theta_cos_model(n, fix_first_theta_hat=True, fix_last_theta_hat=False, theta_delta_max=np.pi/2):

    def con_order(m, i):
        """ Ensures theta_hat[1] < theta_hat[2] < ... < theta_hat[n]."""
        if i == 1:
            return pe.Constraint.Skip
        return m.theta_hat[i] - m.theta_hat[i - 1] >= 1e-6

    def f(theta, theta_hat):
        """Computes the line that is tangent to cos(theta) at theta_hat."""
        return (theta_hat - theta) * pe.sin(theta_hat) + pe.cos(theta_hat)

    def theta_ab(m, i):
        """Computes theta such that f(theta; theta_hat[i-1]) = f(theta; theta_hat[i])."""
        a, b = m.theta_hat[i - 1], m.theta_hat[i]
        return (a * pe.sin(a) + pe.cos(a) - b * pe.sin(b) - pe.cos(b)) / (pe.sin(a) - pe.sin(b))

    def f_theta_ab(m, i):
        """Computes f(theta; theta_hat[i-1]) for theta such that f(theta; theta_hat[i-1]) = f(theta; theta_hat[i])."""
        a, b = m.theta_hat[i - 1], m.theta_hat[i]
        return (pe.sin(a) * pe.cos(b) - pe.cos(a) * pe.sin(b) - (a - b) * pe.sin(a) * pe.sin(b)) / (pe.sin(a) - pe.sin(b))

    def con_error(m, i):
        """Ensures the worst-case error is at least as large as it is at every vertex of the polygonal relaxation."""
        if i == 1:
            return pe.Constraint.Skip
        return m.z >= f_theta_ab(m, i) - pe.cos(theta_ab(m, i))

    m = pe.ConcreteModel()
    m.I_theta_hat = pe.Set(initialize=range(1, n + 1))
    m.z = pe.Var(bounds=(0, 1))
    m.theta_hat = pe.Var(m.I_theta_hat, bounds=(0, theta_delta_max))
    m.obj = pe.Objective(sense=pe.minimize, expr=m.z)
    m.con_order = pe.Constraint(m.I_theta_hat, rule=con_order)
    m.con_error = pe.Constraint(m.I_theta_hat, rule=con_error)
    expr = m.theta_hat[1] == 0 if fix_first_theta_hat else m.z >= f(0, m.theta_hat[1]) - pe.cos(0)
    m.con_first_theta_hat = pe.Constraint(expr=expr)
    expr = m.theta_hat[n] == theta_delta_max if fix_last_theta_hat else m.z >= f(theta_delta_max, m.theta_hat[n]) - pe.cos(theta_delta_max)
    m.con_last_theat_hat = pe.Constraint(expr=expr)
    return m


def get_Theta_cos(n, fix_first_theta_hat=True, fix_last_theta_hat=False, theta_delta_max=np.pi/2):
    m = get_Theta_cos_model(n, fix_first_theta_hat, fix_last_theta_hat, theta_delta_max)
    options = {'first_feasible_solution': True}
    os.environ['NEOS_EMAIL'] = EMAIL
    solver_manager = pe.SolverManagerFactory('neos')
    solver_manager.solve(m, tee=True, solver='knitro')
    return np.array([pe.value(m.theta_hat[i]) for i in m.I_theta_hat])


def obj_breakdown(casestudy, pftype, f_min, f_max, r_hat):

    df_obj = pd.read_csv(os.path.join(CONSOLIDATED_RESULTS_DIR, f'obj-{casestudy}-r{r_hat}.csv'),
                         header=0, index_col=[0, 1, 2, 3])['0']
    df_obj.index.names = ['model', 'casestudy', 'pftype', 'f']
    df_obj = df_obj.loc[(df_obj.index.get_level_values('f') >= f_min) &\
                        (df_obj.index.get_level_values('f') <= f_max)]

    with open(SPECS_FLOOD_YAMLS[casestudy,]) as fh:
        specs = yaml.load(fh, Loader=yaml.Loader)
        probability = pd.Series(specs['probability'])
        probability.index = probability.index.astype(str)

    df_sp = df_obj.loc[(df_obj.index.get_level_values('model') == 'SP') &\
                       (df_obj.index.get_level_values('pftype') == pftype)].copy()
    df_sp.index = df_sp.index.get_level_values('f')

    df_ro = df_obj.loc[(df_obj.index.get_level_values('model') == 'RO') &\
                       (df_obj.index.get_level_values('pftype') == pftype)].copy()
    df_ro.index = df_ro.index.get_level_values('f')

    df_eev = df_obj.loc[(df_obj.index.get_level_values('model') == 'EEV') &\
                        (df_obj.index.get_level_values('pftype') == pftype)].copy()
    df_eev.index = df_eev.index.get_level_values('f')

    df_mmv = df_obj.loc[(df_obj.index.get_level_values('model') == 'MMV') &\
                        (df_obj.index.get_level_values('pftype') == pftype)].copy()
    df_mmv.index = df_mmv.index.get_level_values('f')

    df_ws = df_obj.loc[(df_obj.index.get_level_values('model').str.startswith('WS')) &\
                       (df_obj.index.get_level_values('pftype') == pftype)].to_frame().copy()
    df_ws.columns = ['obj']
    df_ws.reset_index(inplace=True)
    df_ws['omega'] = df_ws['model'].apply(lambda x: x.split('-')[1])
    df_ews = pd.pivot_table(df_ws, index='f', columns='omega', values='obj')\
               .multiply(pd.Series(probability), axis=1).sum(axis=1)
    df_ews.index.name = 'f'
    df_mws = pd.pivot_table(df_ws, index='f', columns='omega', values='obj').max(axis=1)
    df_mws.index.name = 'f'

    df_heur_sp = df_obj.loc[(df_obj.index.get_level_values('model').str.startswith('Heuristic')) &\
                            (df_obj.index.get_level_values('model').str.endswith('(SP)')) &\
                            (df_obj.index.get_level_values('pftype') == pftype)].to_frame().copy()
    df_heur_sp.columns = ['obj']
    df_heur_sp.reset_index(inplace=True)
    df_heur_sp = pd.pivot_table(df_heur_sp, index='f', columns='model', values='obj')
    df_heur_sp.columns = df_heur_sp.columns.map(lambda x: x.strip(' (SP)'))

    df_heur_ro = df_obj.loc[(df_obj.index.get_level_values('model').str.startswith('Heuristic')) &\
                            (df_obj.index.get_level_values('model').str.endswith('(RO)')) &\
                            (df_obj.index.get_level_values('pftype') == pftype)].to_frame().copy()
    df_heur_ro.columns = ['obj']
    df_heur_ro.reset_index(inplace=True)
    df_heur_ro = pd.pivot_table(df_heur_ro, index='f', columns='model', values='obj')
    df_heur_ro.columns = df_heur_ro.columns.map(lambda x: x.strip(' (RO)'))

    return ((df_sp, df_eev, df_ews, df_heur_sp), (df_ro, df_mmv, df_mws, df_heur_ro))


class PyomoSolution:

    def __init__(self, solution):
        self._solution = solution

    def __getitem__(self, variable):
        return self._solution[variable]

    def keys(self):
        return self._solution.keys()

    def to_zip(self, filename, variables=None):
        if variables is None:
            variables = list(self._solution)
        with ZipFile(filename, 'w') as zfh:
            for variable in variables:
                if variable in self._solution:
                    if not type(self._solution[variable]) is pd.Series:
                        zfh.writestr('{}.txt'.format(variable), str(self._solution[variable]))
                    else:
                        zfh.writestr('{}.csv'.format(variable), self._solution[variable].to_csv(header=False))

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
                        solution['objective'] = float(fh.read().strip())
                else:
                    continue
        return cls(solution)

    @classmethod
    def from_solved_instance(cls, instance):
        solution = {}
        for obj in instance.component_objects(pe.Var, active=True):
            if obj.dim() > 1:
                index = pd.MultiIndex.from_tuples(obj.keys())
            else:
                index = list(obj.keys())
            solution[obj.name] = pd.Series(data=[obj[idx].value for idx in index], index=index)
        for obj in instance.component_objects(pe.Objective, active=True):
            solution['objective'] = obj()
        return cls(solution)
