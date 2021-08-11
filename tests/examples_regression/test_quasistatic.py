"""
Created at 21.08.2019
"""

import numpy as np

from conmech.problem_solver import Quasistatic as QuasistaticProblem
from examples.example_quasistatic import QuasistaticSetup

setup = QuasistaticSetup()
runner = QuasistaticProblem(setup, "global optimization")


def test_global_optimization_solver():
    expected_displacement_vector = \
        [-0.11052389326884332, -0.18334981294439903, -0.22349516368943845, -0.24084225308212726, -0.24482530142575815,
         -0.11974107986949582, 0.04034522634248626, -0.08638415150795271, 0.013870068767039412, -0.029274464596135932,
         -0.04436156100281342, -0.03874799756781136, -0.029289873976957495, -0.015496368257876932, 0.0536094994209428,
         -0.14473121207881728, 0.054116689547448846, -0.1377465519422004, 0.05110704481517737, -0.04623677211644456,
         0.15355768963980085, 0.15335220689667983, 0.14728008957040642, 0.12581007130397087, 0.07808333865924907,
         -0.10531586890420432, -0.2589132396690917, -0.45079637936192524, -0.6545354913385428, -0.8560880572010334,
         -0.34726384198889887, -0.16234276818561544, -0.16716806488825947, -0.033508052895117546, -0.03590161969392732,
         -0.652351883981854, -0.4448676868401986, -0.24815498588575652, -0.08413768201257847, -0.7542267476069939,
         -0.7548692986575909, -0.5477173899993708, -0.5495878422108429, -0.3440402789252217, -0.8558889906729761,
         -0.8559953151830919, -0.6520444570141012, -0.44578168186101647, -0.25115414162708916, -0.09501769226549925]
    runner.solving_method = "global optimization"
    results = runner.solve(n_steps=8)
    displacement_vector = results[-1].displacement.T.reshape(1, -1)[0]
    np.testing.assert_array_almost_equal(displacement_vector, expected_displacement_vector, decimal=6)


def test_schur_complement_solver():
    expected_displacement_vector = \
        [-0.07868304304918923, -0.13052853104751239, -0.15910840559867115, -0.17145797330500961, -0.17429354774524258,
         -0.08524484827440978, 0.028722163967603666, -0.061497721938138426, 0.009874239447818442, -0.020840776030175402,
         -0.03158142831321022, -0.027585077403992322, -0.020851746573549122, -0.011032014662176298, 0.03816513666919463,
         -0.1030355728925385, 0.03852620871205764, -0.09806312333296352, 0.03638361052965847, -0.03291640970752091,
         0.10931923668471674, 0.10917295028324045, 0.10485014568237183, 0.08956542371460044, 0.055588291796153884,
         -0.07497539074659129, -0.1843228580635516, -0.3209263440556613, -0.46597021881539125, -0.6094574802419158,
         -0.24722051320781194, -0.11557339196536222, -0.11900856920631045, -0.023854704839400716, -0.025558708368128753,
         -0.46441568567742036, -0.31670564678741553, -0.17666394913914393, -0.059898430273596814, -0.536941403108273,
         -0.5373988434170548, -0.3899253589429035, -0.39125695339579025, -0.24492562580809443, -0.6093157588106195,
         -0.6093914497897699, -0.4641968254073422, -0.31735632900864563, -0.17879907805044676, -0.0676440156673378]
    runner.solving_method = "schur"
    results = runner.solve(n_steps=8)
    displacement_vector = results[-1].displacement.T.reshape(1, -1)[0]
    np.testing.assert_array_almost_equal(displacement_vector, expected_displacement_vector, decimal=6)