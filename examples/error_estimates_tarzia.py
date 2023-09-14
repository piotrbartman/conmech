import numba
import numpy as np
from scipy import interpolate
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import matplotlib.tri as tri
import matplotlib.pylab as pl
from conmech.state.state import TemperatureState

# TODO #99


def compare(ref: TemperatureState, sol: TemperatureState):
    tt = 0
    x = sol.body.mesh.initial_nodes[:, 0]
    y = sol.body.mesh.initial_nodes[:, 1]

    soltri = tri.Triangulation(x, y, triangles=sol.body.mesh.elements)
    thi = tri.LinearTriInterpolator(soltri, sol.temperature)

    for element in ref.body.mesh.elements:
        x0 = ref.body.mesh.initial_nodes[element[0]]
        x1 = ref.body.mesh.initial_nodes[element[1]]
        x2 = ref.body.mesh.initial_nodes[element[2]]
        t0 = ref.temperature[element[0]]
        t1 = ref.temperature[element[1]]
        t2 = ref.temperature[element[2]]
        t = t0 + t1 + t2

        tdx, tdy = calculate_dx_dy(x0, t0, x1, t1, x2, t2)
        th0 = np.ma.getdata(thi(*x0))
        th1 = np.ma.getdata(thi(*x1))
        th3 = np.ma.getdata(thi(*x2))
        thdx, thdy = calculate_dx_dy(x0, th0, x1, th1, x2, th3)

        th = thi(*x0) + thi(*x1) + thi(*x2)
        tt += ((t - th) ** 2 + (tdx - thdx) ** 2 + (tdy - thdy) ** 2) ** 0.5
    return tt


@numba.njit()
def calculate_dx_dy(x0, u0, x1, u1, x2, u2):
    a1 = x1[0] - x0[0]
    b1 = x1[1] - x0[1]
    c1 = u1 - u0
    a2 = x2[0] - x0[0]
    b2 = x2[1] - x0[1]
    c2 = u2 - u0
    a = b1 * c2 - b2 * c1
    b = a2 * c1 - a1 * c2
    c = a1 * b2 - b1 * a2
    dx = a / c
    dy = b / c
    return dx, dy
