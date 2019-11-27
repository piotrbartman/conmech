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


def integral_phi_i_x_1__phi_j_x_1(w: np.ndarray, edges, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2
    if etype == UP_TRIANGLE:
        type_side(edges, element, w, efield)

    if etype == RIGHT_TRIANGLE:
        type_up(edges, element, w, efield)

    if etype == BOTTOM_TRIANGLE:
        type_side(edges, element, w, efield)

    if etype == LEFT_TRIANGLE:
        type_up(edges, element, w, efield)


def integral_phi_i_x_2__phi_j_x_2(w: np.ndarray, edges, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2

    if etype == UP_TRIANGLE:
        type_up(edges, element, w, efield)

    if etype == RIGHT_TRIANGLE:
        type_side(edges, element, w, efield)

    if etype == LEFT_TRIANGLE:
        type_side(edges, element, w, efield)

    if etype == BOTTOM_TRIANGLE:
        type_up(edges, element, w, efield)


def integral_phi_i_x_1__phi_j_x_2(w: np.ndarray, edges, element, etype):

    efield = 1/4 * edges.length[element[EDGE_0]] ** 2

    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]
    c = edges[element[EDGE_1]][0]
    value = (1 / edges.length[element[EDGE_0]]) ** 2 * efield

    if etype == UP_TRIANGLE:
        w[a, b] += - value
        w[b, a] += value
        w[a, c] += 2 * value
        w[b, c] += - 2 * value
        w[a, a] += - value
        w[b, b] += value

    if etype == RIGHT_TRIANGLE:
        w[a, b] += value
        w[b, a] += - value
        w[c, a] += 2 * value
        w[c, b] += - 2 * value
        w[a, a] += - value
        w[b, b] += value

    if etype == BOTTOM_TRIANGLE:
        w[a, b] += value
        w[b, a] += - value
        w[a, c] += - 2 * value
        w[b, c] += 2 * value
        w[a, a] += value
        w[b, b] += - value

    if etype == LEFT_TRIANGLE:
        w[a, b] += - value
        w[b, a] += value
        w[c, a] += - 2 * value
        w[c, b] += 2 * value
        w[a, a] += value
        w[b, b] += - value


def type_up(edges, element, w, efield):
    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]

    w_a_b = - (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    w[a, b] += w_a_b
    w[b, a] += w_a_b

    w_a_a = (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    w[a, a] += w_a_a
    w[b, b] += w_a_a


def type_side(edges, element, w, efield):
    a = edges[element[EDGE_0]][0]
    b = edges[element[EDGE_0]][1]
    c = edges[element[EDGE_1]][0]

    w_a_b = (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    w[a, b] += w_a_b
    w[b, a] += w_a_b

    w_a_a = w_a_b
    w[a, a] += w_a_a
    w[b, b] += w_a_a

    w[c, c] += (-2 / edges.length[element[EDGE_0]]) ** 2 * efield

    w_c_b = (-2) * (1 / edges.length[element[EDGE_0]]) ** 2 * efield
    w[c, b] += w_c_b
    w[b, c] += w_c_b

    w_c_a = w_c_b
    w[c, a] += w_c_a
    w[a, c] += w_c_a


def w11(mesh: Mesh):
    w11 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_1__phi_j_x_1(w11, mesh.edge, element, mesh.element.type[i])
    return w11


def w22(mesh: Mesh):
    w22 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_2__phi_j_x_2(w22, mesh.edge, element, mesh.element.type[i])
    return w22


def w12(mesh: Mesh):
    w12 = np.zeros((len(mesh.point), len(mesh.point)))
    for i, element in enumerate(mesh.element):
        integral_phi_i_x_1__phi_j_x_2(w12, mesh.edge, element, mesh.element.type[i])
    return w12

