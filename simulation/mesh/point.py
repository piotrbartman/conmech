"""
Created at 22.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class Point:
    LEFT_BOTTOM_CORNER = 0
    LEFT_SIDE = 1
    LEFT_TOP_CORNER = 2
    TOP = 3
    RIGHT_TOP_CORNER = 4
    RIGHT_SIDE = 5
    RIGHT_BOTTOM_CORNER = 6
    BOTTOM = 7
    NORMAL_MIDDLE = 8
    CROSS = 9

    @staticmethod
    def gradients(point):
        dx = np.array([1., 1., -1., -1., -1., -1., 1., 1.]) * 0.5  # normal dx
        dy = np.array([-1., -1., -1., -1., 1., 1., 1., 1.]) * 0.5
        
        if point[2] == Point.TOP:
            f = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        elif point[2] == Point.RIGHT_TOP_CORNER:
            f = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        elif point[2] == Point.RIGHT_SIDE:
            f = np.array([1, 1, 0, 0, 0, 0, 1, 1])
        elif point[2] == Point.RIGHT_BOTTOM_CORNER:
            f = np.array([1, 1, 0, 0, 0, 0, 0, 0])
        elif point[2] == Point.BOTTOM:
            f = np.array([1, 1, 1, 1, 0, 0, 0, 0])
        elif point[2] == Point.NORMAL_MIDDLE:
            f = np.array([1, 1, 1, 1, 1, 1, 1, 1])
        elif point[2] == Point.CROSS:
            f = np.array([1, 1, 1, 1, 0, 0, 0, 0])  # only 4 used
            dx = np.array([1., 0., -1., 0., 0., 0., 0., 0.])  # cross dx
            dy = np.array([0., -1., 0., 1., 0., 0., 0., 0.])
        else:
            raise ValueError

        result = (f * dx, f * dy)
        return result

