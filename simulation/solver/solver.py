"""
Created at 21.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
import scipy.optimize
import numba


class Solver:

    def __init__(self, mesh, alpha, regular_dphi, b, rho):

        self.mesh = mesh

        self.alpha = alpha
        self.regular_dphi = regular_dphi
        self.b = b
        self.rho = rho

        self.u = None
        self.ub = None

        self.B = None
        self.F = None

        # TODO move to the setup
        self.precision_coefficient = 1e8

    def solve(self, start_vector=None, verbose=True, quality=1):
        u_vector_shape = (self.mesh.ind_num,)
        if start_vector is None:
            u_vector = np.zeros(u_vector_shape)
        else:
            assert start_vector.shape == u_vector_shape
            u_vector = start_vector

        while True:
            u_vector = scipy.optimize.fsolve(self.f, u_vector)
            quality_inv = np.linalg.norm(self.f(u_vector))
            if quality_inv < quality ** -1:
                if verbose:
                    if quality_inv == 0:
                        print("Found exact solution.")
                    else:
                        print(f"Quality = {quality_inv ** -1} > {quality} is acceptable.")
                break
            else:
                if verbose:
                    print(f"Quality = {quality_inv ** -1} is too low, trying again...")

    # TODO cleanup
    @staticmethod
    @numba.jit()
    def numba_JZu(ind_num, Edges, Points, u, BorderEdgesD, BorderEdgesN, BorderEdgesC, regular_dphi, b, rho):
        result = np.zeros(ind_num)

        for i in range(ind_num):
            for e in range(-BorderEdgesD - BorderEdgesN - BorderEdgesC,
                           -BorderEdgesD - BorderEdgesN):
                e1 = Edges[e][0]
                e2 = Edges[e][1]
                if i == e1 or i == e2:
                    umL = _u_at_middle(e1, e2, ind_num, u)

                    p1 = Points[e1][:2]
                    p2 = Points[e2][:2]

                    L = _distance(p1, p2)
                    result[i] += L * 0.5 * regular_dphi(umL, b, rho)

        return result

    def JZu(self):
        JZu = Solver.numba_JZu(self.mesh.ind_num, self.mesh.Edges, self.mesh.Points, self.u,
                               self.mesh.borders["Dirichlet"], self.mesh.borders["Neumann"], self.mesh.borders["Contact"],
                               self.regular_dphi, self.b, self.rho)
        return JZu

    # TODO cleanup
    @staticmethod
    @numba.jit()
    def numba_Bu1(B, u):
        result = np.dot(B, u)
        return result

    def Bu1(self):
        result = Solver.numba_Bu1(self.B, self.u)
        return result

    def condition(self):
        return self.alpha * self.JZu()

    def f(self, u_vector):
        self.u = u_vector

        X = self.Bu1() \
            + self.condition() \
            - self.F.F

        return self.precision_coefficient * X


@numba.njit()
def _distance(p1, p2):
    return np.sqrt(np.sum((p1 - p2) ** 2))


@numba.njit()
def _u_at_middle(e1, e2, ind_num, u):
    result = 0
    if e1 < ind_num:
        result += u[e1] * 0.5
    if e2 < ind_num:
        result += u[e2] * 0.5
    return result
