"""
Created at 30.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class Point:
    def __init__(self, dimension, n):
        self.coordinates = np.empty((n, dimension))

    def __getitem__(self, item):
        return self.coordinates[item]

    def __len__(self):
        return len(self.coordinates)
