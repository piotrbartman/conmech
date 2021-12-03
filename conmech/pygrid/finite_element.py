
from conmech.pygrid import point


class FiniteElement:

    def __init__(self, point_1: "point.Point", point_2: "point.Point", point_3: "point.Point", /):
        self.points: tuple = (point_1, point_2, point_3)
        self.area = .5 * abs((point_2.x - point_1.x) * (point_3.y - point_1.y)
                           - (point_2.y - point_1.y) * (point_3.x - point_1.x))
