"""
Created at 13.11.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np
from new.constants import *
from new.mesh import Mesh

"""
 a    b
b  ---  b
 |\   /|
 | \ / |
 |  c  |
 | / \ |
 |/   \|
a  ---  a
 a     b
"""


def integral_phi_i_x_1__phi_j_x_1(w: np.ndarray, edges, points, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2
    if etype == UP_TRIANGLE:
        type_side(edges, points, element, w, efield)

    if etype == RIGHT_TRIANGLE:
        type_up(edges, points, element, w, efield)

    if etype == BOTTOM_TRIANGLE:
        type_side(edges, points, element, w, efield)

    if etype == LEFT_TRIANGLE:
        type_up(edges, points, element, w, efield)


def integral_phi_i_x_2__phi_j_x_2(w: np.ndarray, edges, points, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2

    if etype == UP_TRIANGLE:
        type_up(edges, points, element, w, efield)

    if etype == RIGHT_TRIANGLE:
        type_side(edges, points, element, w, efield)

    if etype == LEFT_TRIANGLE:
        type_side(edges, points, element, w, efield)

    if etype == BOTTOM_TRIANGLE:
        type_up(edges, points, element, w, efield)


def integral_phi_i_x_1__phi_j_x_2(w: np.ndarray, edges, points, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2

    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]
    c = edges[element[EDGE_1]][0]
    a_ok = not points.dirichlet[a]
    b_ok = not points.dirichlet[b]
    c_ok = not points.dirichlet[c]
    value = (1 / edges.length[element[EDGE_0]]) ** 2 * efield

    if etype == UP_TRIANGLE:
        if a_ok and b_ok: w[a, b] += - value
        if a_ok and b_ok: w[b, a] += value
        if a_ok and c_ok: w[a, c] += 2 * value
        if c_ok and b_ok: w[b, c] += - 2 * value
        if a_ok: w[a, a] += - value
        if b_ok: w[b, b] += value

    if etype == RIGHT_TRIANGLE:
        if a_ok and b_ok: w[a, b] += value
        if a_ok and b_ok: w[b, a] += - value
        if a_ok and c_ok: w[c, a] += 2 * value
        if c_ok and b_ok: w[c, b] += - 2 * value
        if a_ok: w[a, a] += - value
        if b_ok: w[b, b] += value

    if etype == BOTTOM_TRIANGLE:
        if a_ok and b_ok: w[a, b] += value
        if a_ok and b_ok: w[b, a] += - value
        if a_ok and c_ok: w[a, c] += - 2 * value
        if c_ok and b_ok: w[b, c] += 2 * value
        if a_ok: w[a, a] += value
        if b_ok: w[b, b] += - value

    if etype == LEFT_TRIANGLE:
        if a_ok and b_ok: w[a, b] += - value
        if a_ok and b_ok: w[b, a] += value
        if a_ok and c_ok: w[c, a] += - 2 * value
        if c_ok and b_ok: w[c, b] += 2 * value
        if a_ok: w[a, a] += value
        if b_ok: w[b, b] += - value


def type_up(edges, points, element, w, efield):
    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]
    a_ok = not points.dirichlet[a]
    b_ok = not points.dirichlet[b]

    if a_ok and b_ok:
        w_a_b = - (1 / edges.length[element[EDGE_0]]) ** 2 * efield
        w[a, b] += w_a_b
        w[b, a] += w_a_b

    w_a_a = (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    if a_ok: w[a, a] += w_a_a
    if b_ok: w[b, b] += w_a_a


def type_side(edges, points, element, w, efield):
    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]
    c = edges[element[EDGE_1]][0]
    a_ok = not points.dirichlet[a]
    b_ok = not points.dirichlet[b]
    c_ok = not points.dirichlet[c]

    w_a_b = (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    if a_ok and b_ok:
        w[a, b] += w_a_b
        w[b, a] += w_a_b

    w_a_a = w_a_b
    if a_ok: w[a, a] += w_a_a
    if b_ok: w[b, b] += w_a_a
    if c_ok: w[c, c] += (-2 / edges.length[element[EDGE_0]]) ** 2 * efield

    if c_ok:
        w_c_b = (-2) * (1 / edges.length[element[EDGE_0]]) ** 2 * efield
        if b_ok:
            w[c, b] += w_c_b
            w[b, c] += w_c_b

        if a_ok:
            w_c_a = w_c_b
            w[c, a] += w_c_a
            w[a, c] += w_c_a


def w11(mesh: Mesh):
    w11 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_1__phi_j_x_1(w11, mesh.edge, mesh.point, element, mesh.element.type[i])
    return w11


def w22(mesh: Mesh):
    w22 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_2__phi_j_x_2(w22, mesh.edge, mesh.point, element, mesh.element.type[i])
    return w22


def w12(mesh: Mesh):
    w12 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_1__phi_j_x_2(w12, mesh.edge, mesh.point, element, mesh.element.type[i])
    return w12

