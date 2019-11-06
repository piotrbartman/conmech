"""
Created at 30.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class Edge:
    def __init__(self, dimension, n):
        self.points = np.empty((n, dimension))

    def __getitem__(self, item):
        return self.points[item]

    def __len__(self):
        return len(self.points)
