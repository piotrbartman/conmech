"""
Created at 29.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

from new.point import Point
from new.edge import Edge
from new.element import Element
from new.smart_array import SmartArray
import numpy as np


class Mesh:

    # TODO: rename point -> points etc.
    def __init__(self, point: Point, edge: Edge, element: Element,  subarea: dict):
        self.point = point
        self.edge = edge
        self.element = element
        self.subarea = subarea
        # self.independent = SmartArray(self.point, idx)
