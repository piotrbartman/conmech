"""
Created at 29.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

from new.point import Point
from new.edge import Edge
from new.element import Element


class Mesh:
    def __init__(self, point: Point, edge: Edge, element: Element,  subarea: dict):
        self.point = point
        self.edge = edge
        self.element = element
        self.subarea = subarea
