"""
Created at 30.10.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""

import numpy as np


class Element:
    def __init__(self, dimension, n):
        # TODO: dict?
        self.edges = np.empty((n, dimension + 1), dtype=np.int32)
        self.type = np.empty(n, dtype=np.int32)
        self.fields = np.empty(n)

    def __getitem__(self, item):
        return self.edges[item]

    def __len__(self):
        return len(self.edges)
