# coding=windows-1251

"""
Module for symbolic derivation of equations of motion of mechanical systems.
"""

from ctypes import ArgumentError
from symbol import except_clause
from sympy.core.numbers import Zero
from sympy import (Symbol, Function, symbols, sin, cos, tan, cot,
    solve, zeros, Derivative as D, diff, Eq, collect, Matrix, pprint,
    simplify, radsimp, fraction, together, expand)
from sympy.printing.str import StrPrinter

t = Symbol('t')

class SymPrinter(StrPrinter):

    def _print_Function(self, e):
        """
        Print ui(t) as ui, where is i is an index number of the generalized
        speed.
        """
        if hasattr(e, 'is_gc'):
            return str(e.func)
        else:
            return StrPrinter().doprint(e)

    def _print_Symbol(self, e):
        """
        Print ui(t) as ui, where is i is an index number of the generalized
        speed.
        """
        if hasattr(e, 'is_gc'):
            return str(e.func)
        else:
            return StrPrinter().doprint(e)

def printm(expr):
    op = lambda expr: str(expr).replace('**', '^').replace('D(', 'diff(')
    if isinstance(expr, list):
        for e in expr:
            print op(e)
    else:
        print op(expr)

def GeneralizedCoordinate(name, constant=False):
    gc = Symbol(name)(Symbol('t'))
    gc.is_gc = True
    if constant:
        gc.fdiff = lambda argindex: 0
    gc.__repr__ = lambda self: SymPrinter().doprint(self)
    gc.__str__ = lambda self: SymPrinter().doprint(self)
    return gc

def gcs(s, number=1, from_num=1):
    gc_list = [GeneralizedCoordinate(s[0] + str(i)) for i in range(from_num, from_num + number)]
    gcd_list = [gc.diff(t) for gc in gc_list]
    gc2d_list = [gcd.diff(t) for gcd in gcd_list]
    if number == 1:
        return (gc_list[0], gcd_list[0], gc2d_list[0])
    else:
        return (gc_list, gcd_list, gc2d_list)

def pdiff(expr, arg):
    """
    Partial derivative.
    """
    s = Symbol('_s_')
    op = lambda f, x: diff(f.subs(x, s), s).subs(s, x)
    if isinstance(arg, list):
        return [op(expr, i) for i in arg]
    else:
        return op(expr, arg)

def lagrange_equations_lhs(lagrangian, q):
    """
    Left-hand side of Lagrange equation.
    """
    op = lambda f, x: diff(pdiff(f, diff(x, t)), t) - pdiff(f, x)
    if isinstance(q, list):
        return [op(lagrangian, qi) for qi in q]
    elif isinstance(q, Matrix):
        res = zeros(q.shape)
        for i in range(res.shape[0]):
            for j in range(res.shape[1]):
                res[i, j] = op(lagrangian, q[i, j])
        return res
    else:
        return op(lagrangian, q)

def solve_slae_by_gauss(eqns, vars):
    """
    Solves SLAE by Gauss's method.
    """
    n = len(vars)
    res = list(vars)
    for i in range(n):
        tmp = solve(eqns[i], vars[i])
        if len(tmp):
            res[i] = tmp[0]
        for j in range(i+1, n):
            eqns[j] = eqns[j].subs(vars[i], res[i])
    for i in range(n-1, -1, -1):
        for j in range(i-1, -1, -1):
            res[j] = res[j].subs(vars[i], res[i])
    return dict(zip(vars, res))

def normalize(equations, expanded=True):
    """
    Normalization of differential equations.
    """
    if not isinstance(equations, dict):
        raise ArgumentError, "equations must be dict"
    coeffs = {}
    neqns = {}
    vars = []
    i = 1
    for k in equations.keys():
        tmp = equations[k]
        coeffs[k] = {}
        neqns[k] = 0
        j = 1
        for a in equations.keys():
            c = pdiff(tmp, a)
            if not isinstance(c, Zero) and not c is 0:
                c_name = Symbol('C' + str(i) + str(j))
                coeffs[k][c_name] = c
                tmp -= c*a
                neqns[k] += c_name * a
                vars.append(a)
                j += 1
        if tmp is not 0:
            c_name = Symbol('C' + str(i) + '0')
            coeffs[k][c_name] = tmp
            neqns[k] += c_name
        i += 1
    neqns = solve_slae_by_gauss(neqns.values(), vars)
    for k in neqns:
        for ck in coeffs[k]:
            if expanded:
                neqns[k] = expand(neqns[k].subs(ck, coeffs[k][ck]))
            else:
                neqns[k] = neqns[k].subs(ck, coeffs[k][ck])
    return neqns

class MechanicalFrame:

    def __init__(self, name=''):
        """
        Initialization
        """
        self.name = name

        # Joint coordinates
        self.q_dim = 0
        self.q_list = []
        self.u_list = []
        self.a_list = []
        self.u_dependent = []
        self.u_independent = []
        self.u_names_dict = {}
        self.a_names_dict = {}

        # Parameters and forces
        self.parameter_list = []
        self.joint_forces = {}
        self.control = {}

        # Equilibrium point
        self.q0_list = []
        self.u0_list = []

        self.hc_eqns = []
        self.dhc_eqns = []
        self.dhc_matrix = Matrix()
        #self.template = MatrixTemplate()

    def add_coordinates(self, string, number=1):
        """
        Declares generalized coordinates, velocities and accelerations.
        """
        from_num = self.q_dim + 1 if self.q_dim > 0 else 1
        q_list, qdot_list, q2dot_list = gcs(string, number, from_num)

        # Coordinates
        try:
            for q in q_list:
                self.q_list.append(q)
        except TypeError:
            self.q_list.append(q_list)
        
        # Speeds
        try:
            for qd in qdot_list:
                self.u_list.append(qd)
                self.u_independent.append(qd)
        except TypeError:
            if qdot_list:
                self.u_list.append(qdot_list)
            else:
                self.u_list.append(q_list.diff(t))

        self.u_names_dict = dict(zip(self.u_list, [Symbol('qd' + str(i+1)) for i in range(len(self.q_list))]))

        # Accelerations
        try:
            for q2d in q2dot_list:
                self.a_list.append(q2d)
        except TypeError:
            if q2dot_list:
                self.a_list.append(q2dot_list)
            else:
                self.a_list.append(qdot_list.diff(t))

        self.a_names_dict = dict(zip(self.a_list, [Symbol('q2d' + str(i+1)) for i in range(len(self.q_list))]))

        self.q_dim += number

        return q_list, qdot_list, q2dot_list

    def add_parameters(self, string):
        """
        Declares parameters of mechanical system.
        """
        params = list(symbols(string))
        for p in params:
            self.parameter_list.append(p)
        return params

    def add_joint_forces(self, joint_forces):
        """
        Declares joint forces of mechanical systems.
        """
        for q, f in joint_forces.iteritems():
            self.joint_forces[q] = f

    def add_control_forces(self, control_forces):
        """
        Declares control forces.
        """
        for q, f in control_forces.iteritems():
            self.control_forces[q] = f

    def form_constraints_matrix(self, dhc_eqns, u_dep_list):
        """
        Calculates matrix of differential constraints reference to self.u_list
        or to u_list if it's presented.
        """
        self.dhc_eqns = dhc_eqns
        self.u_dependent = u_dep_list
        self.u_independent = [u for u in self.u_list if u not in u_dep_list]
        m = len(self.dhc_eqns)
        k = len(self.u_independent)
        self.dhc_matrix = zeros([m, k])
        for i in range(m):
            tmp = self.dhc_eqns[i]
            for j in range(k):
                self.dhc_matrix[i, j] = pdiff(tmp, self.u_independent[j])
                tmp -= self.dhc_matrix[i, j] * self.u_independent[j]
        return self.dhc_matrix

    def form_lagranges_equations(self):
        """
        Calculates Lagrange's equations. The method uses undetermined
        multipliers if the system has differential constraints.
        """
        m = len(self.q_list)
        n = len(self.dhc_eqns)

        # Lagrange's multipliers
        self.lambda_list = [Symbol('lambda' + str(i)) for i in range(n)]

        # Equations
        eqns = {}
        for q in self.q_list:
            tmp = lagrange_equations_lhs(self.lagrangian, q) - self.joint_forces.get(q, 0)
            for j in range(n):
                tmp -= self.lambda_list[j] * self.dhc_matrix[j, i]
            eqns[diff(self.q_list[i], t, t)] = tmp

        # Excluding Lagrange's multipliers
        if len(self.lambda_list):
            u_dep = solve(self.dhc_eqns, self.u_dependent);
            a_dep = [diff(u, t) for u in u_dep]
            for k in eqns.keys():
                eqns[k] = eqns[k].subs(zip(self.u_list, u_dep))
                eqns[k] = eqns[k].subs(zip(self.a_list, a_dep))
                eqns[k] = simplify(eqns[k])

        self.lagrange_eqnuations = eqns
        return self.lagrange_eqnuations

    def form_shulgins_equations(self, normalized=False, packed=False, expanded=True):
        """
        Calculates Shulgin's equations for systems with (or without)
        redundant coordinates.
        """
        m = len(self.u_independent)
        n = len(self.u_dependent)
        assert (n, m) == self.dhc_matrix.shape

        # Let's calculate reduced lagrangian
        L = self.lagrangian.subs(zip(self.u_dependent, self.dhc_eqns))
        self.reduced_lagrangian = L
        
        q_indep = [q for q in self.q_list if q.diff(t) in self.u_independent]
        q_dep = [q for q in self.q_list if q.diff(t) in self.u_dependent]

        # Calculating equations
        eqns = {}
        i = 0
        for q in q_indep:
            tmp = lagrange_equations_lhs(L, q) - self.joint_forces.get(q, 0)
            j = 0
            for q1 in q_dep:
                tmp -= self.dhc_matrix[j, i] * pdiff(L, q1)
                j += 1
            eqns[diff(q, t, t)] = tmp
            i += 1

        if normalized:
            eqns = normalize(eqns, expanded)
            
        # Adding diff. holonomic constraints equations
        for i in range(n):
            if normalized:
                eqns[self.u_dependent[i]] = self.dhc_eqns[i]
            else:
                eqns[self.u_dependent[i]] = self.u_dependent[i] - self.dhc_eqns[i]

        self.shulgins_equations = eqns
        return self.shulgins_equations

    def define_equilibrium_point(self, motion_equations={}):
        """
        Defines equilibrium symbols and determines positional and cyclic
        coordinates if paramater motion_equations is passed.
        """
        self.q0 = dict([(q, 0) for q in self.q_list])
        self.u0 = dict([(u, 0) for u in self.u_list])
        self.a0 = dict([(a, 0) for a in self.a_list])

        if not isinstance(motion_equations, dict):
            raise ArgumentError, "motion_equations must be dict"

        for i in range(self.q_dim):
            u = self.u_list[i]
            q = self.q_list[i]
            u0 = Symbol('u0' + str(i+1)) if u in self.u_independent else 0
            self.u0[u] = u0
            self.q0[q] = u0*t + Symbol('q0' + str(i+1))
            # Check for positional coordinates
            if len(motion_equations):
                for k in motion_equations.keys():
                    tmp = pdiff(motion_equations[k], q)
                    if not tmp is 0 and not isinstance(tmp, Zero):
                        # Looks like positional coordinate
                        self.q0[q] = self.q0[q].subs(u0, 0)
                        self.u0[u] = 0
                        break
        return self.q0, self.u0

    def form_equilibrium_manifold_equations(self, motion_equations={}, point={}):
        """
        Calculates equations of equilibrium.
        """
        self.manifold_equations = {}
        if len(motion_equations):
            if not isinstance(motion_equations, dict):
                raise ArgumentError, "motion_equations must be dict"
            for k in motion_equations.keys():
                tmp = motion_equations[k]
                try:
                    tmp = tmp.subs(self.a0).subs(self.u0).subs(self.q0)
                    if len(point):
                        tmp = simplify(tmp.subs(point))
                finally:
                    self.manifold_equations[k] = tmp
        return self.manifold_equations

    def form_perturbed_equations(self, motion_equations={}, manifold={}):
        """
        Calculates equations of perturbed motion.
        """
        self.perturbed_equations = {}
        if not isinstance(motion_equations, dict):
            raise ArgumentError, "motion_equations must be dict"
                
        # Lets create array of perturbed variables
        self.x = {}
        self.x_list = []
        for i in range(self.q_dim):
            x_i = Symbol('x' + str(i + 1))(t)
            self.x_list.append(x_i)
            self.x[self.q_list[i]] = x_i
        i = 1
        for k in motion_equations.keys():
            if k in self.a_list:
                x_i = Symbol('x' + str(self.q_dim + i))(t)
                self.x_list.append(x_i)
                self.x[self.find_velocity(k)] = x_i
                i += 1
        
        for k in motion_equations.keys():
            tmp = motion_equations[k]
            for q in self.q_list:
                u = q.diff(t)
                if u in self.x:
                    tmp = tmp.subs(u, self.x.get(u) + self.u0.get(u))
                tmp = tmp.subs(q, self.x.get(q) + self.q0.get(q))
            # dx/dt = F(q0 + x) - F(q0)
            self.perturbed_equations[k] = tmp - manifold.get(k, 0)

        return self.perturbed_equations

    def form_first_approximation_equations(self, motion_equations={}, q0={}, \
        u0={}, params={}, simplified=False):
        """
        Calculates first approxamation equitions.
        """
        self.fa_equations = {}
        if not isinstance(motion_equations, dict):
            raise ArgumentError, "motion_equations must be dict"

        n = len(motion_equations)
        m = len(self.u_dependent)

        # Point x = 0
        x_0 = [(x, 0) for x in self.x.values()]

        # OK, now lets find them
        for k in motion_equations.keys():
            tmp = 0
            for q, x in self.x.iteritems():
                tmp += pdiff(motion_equations[k], x).subs(x_0) * x
            if len(q0):
                tmp = tmp.subs(q0)
            if len(u0):
                tmp = tmp.subs(u0)
            if len(params):
                tmp = tmp.subs(params)
            if simplified:
                tmp = simplify(tmp)
            x1 = self.x.get(self.find_coordinate(k))
            if not k in self.a_list:
                self.fa_equations[x1.diff(t)] = tmp
            else:
                x2 = self.x.get(self.find_velocity(k), 0)
                self.fa_equations[x1.diff(t)] = x2
                self.fa_equations[x2.diff(t)] = tmp
            
        return self.fa_equations

    def create_matrix_of_coeff(self, expr_list, vars):
        n = len(expr_list)
        m = len(vars)
        matrix = zeros([n, m])
        for i in range(0, n):
            for j in range(0, m):
                matrix[i, j] = pdiff(expr_list[i], vars[j])
        return matrix

    def find_coordinate(self, var):
        for q in self.q_list:
            if diff(q, t) == var or diff(diff(q, t), t) == var:
                return q

    def find_velocity(self, var):
        for u in self.u_list:
            if diff(u, t) == var:
                return u

    def set_dhc_matrix(self, matrix, u_dep_list):
        """
        Declares differentiated holonomic constraints matrix and lists of
        dependent velocities.
        """
        self.u_dependent = u_dep_list
        self.u_independent = [u for u in self.u_list if u not in u_dep_list]
        self.dhc_matrix = matrix

    def set_lagrangian(self, expr):
        self.lagrangian = expr

    def set_hc_eqns(self, hc_eqns_list):
        """
        Sets holonomic constraints equations.
        """
        self.hc_eqns = hc_eqns_list

    def set_dhc_eqns(self, dhc_eqns_list):
        """
        Sets differentiated holonomic constraints equations
        """
        self.dhc_eqns = dhc_eqns_list

