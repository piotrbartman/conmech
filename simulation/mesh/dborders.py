"""
Created at 03.09.2019

@author: Piotr Bartman
"""


class DBorders:
    DIRICHLET = "Dirichlet"
    NEUMANN = "Neumann"
    CONTACT = "Contact"

    def __init__(self, left, top, right, bottom):
        self.left = left
        self.top = top
        self.right = right
        self.bottom = bottom
