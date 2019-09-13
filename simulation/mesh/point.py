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
    def is_inside(point):
        return point[2] == Point.CROSS or point[2] == Point.NORMAL_MIDDLE

    @staticmethod
    def gradients(point):
        if point[2] == Point.CROSS:
            dx = np.array([1., 0., -1., 0., 0., 0., 0., 0.])
            dy = np.array([0., -1., 0., 1., 0., 0., 0., 0.])
        else:   # normal dx
            dx = np.array([1., 1., -1., -1., -1., -1., 1., 1.]) * 0.5
            dy = np.array([-1., -1., -1., -1., 1., 1., 1., 1.]) * 0.5
        
        f = Point.triangles(point)

        result = (f * dx, f * dy)
        return result

    @staticmethod
    def triangles(point):
        if point[2] == Point.LEFT_BOTTOM_CORNER:
            result = np.array([0, 0, 1, 1, 0, 0, 0, 0])
        elif point[2] == Point.LEFT_SIDE:
            result = np.array([0, 0, 1, 1, 1, 1, 0, 0])
        elif point[2] == Point.LEFT_TOP_CORNER:
            result = np.array([0, 0, 0, 0, 1, 1, 0, 0])
        elif point[2] == Point.TOP:
            result = np.array([0, 0, 0, 0, 1, 1, 1, 1])
        elif point[2] == Point.RIGHT_TOP_CORNER:
            result = np.array([0, 0, 0, 0, 0, 0, 1, 1])
        elif point[2] == Point.RIGHT_SIDE:
            result = np.array([1, 1, 0, 0, 0, 0, 1, 1])
        elif point[2] == Point.RIGHT_BOTTOM_CORNER:
            result = np.array([1, 1, 0, 0, 0, 0, 0, 0])
        elif point[2] == Point.BOTTOM:
            result = np.array([1, 1, 1, 1, 0, 0, 0, 0])
        elif point[2] == Point.NORMAL_MIDDLE:
            result = np.array([1, 1, 1, 1, 1, 1, 1, 1])
        elif point[2] == Point.CROSS:
            result = np.array([1, 1, 1, 1, 0, 0, 0, 0])  # only 4 used
        else:
            raise ValueError

        return result

