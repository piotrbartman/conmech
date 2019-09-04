"""
Created at 03.09.2019

@author: Piotr Bartman
"""


class DCell:
    def __init__(self, x_height, y_height):
        self.x_height = x_height
        self.y_height = y_height
        self.hypotenuse_length = (x_height ** 2 + y_height ** 2) ** 0.5
        self.area = x_height * y_height / 2
