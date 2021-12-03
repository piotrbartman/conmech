
from conmech.pygrid import point


class Edge:

    def __init__(self, start_point: "point.Point", stop_point: "point.Point"):
        self.start_point = start_point
        self.stop_point = stop_point

    def length(self):
        return ((self.start_point.x - self.stop_point.x) ** 2 + (self.start_point.y - self.stop_point.y) ** 2) ** .5
