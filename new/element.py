"""
Created at 30.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class Element:
    def __init__(self, dimension, n):
        self.edges = np.empty((n, dimension + 1))

    def __getitem__(self, item):
        return self.edges[item]

    def __len__(self):
        return len(self.edges)
