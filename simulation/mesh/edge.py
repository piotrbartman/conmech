"""
Created at 22.08.2019

@author: Michał Jureczka
@author: Piotr Bartman
"""


class Edge:
    # TODO:  documentation
    # i, j, type: (always i<j on plane)
    NON_EXISTENT = 0  # no edge
    HORIZONTAL = 1  # from normal go right to normal
    VERTICAL = 2  # from normal go up to normal,
    TOP = 3  # from normal go right and up to cross
    X_TOP = 4  # from cross go right and up to normal,
    BOTTOM = 5  # from normal go right and down to cross
    X_BOTTOM = 6  # from cross go right and down to normal

    @staticmethod
    def c(edge):
        if edge[2] == Edge.HORIZONTAL:
            c1i = 3
            c1j = 0
            c2i = 4
            c2j = 7
        elif edge[2] == Edge.VERTICAL:
            c1i = 1
            c1j = 6
            c2i = 2
            c2j = 5
        elif edge[2] == Edge.TOP:
            c1i = 2
            c1j = 0
            c2i = 3
            c2j = 3
        elif edge[2] == Edge.X_TOP:
            c1i = 1
            c1j = 7
            c2i = 2
            c2j = 6
        elif edge[2] == Edge.BOTTOM:
            c1i = 4
            c1j = 1
            c2i = 5
            c2j = 0
        elif edge[2] == Edge.X_BOTTOM:
            c1i = 2
            c1j = 1
            c2i = 3
            c2j = 0
        else:
            print("ups")
            c1i = -1
            c1j = -1
            c2i = -1
            c2j = -1

        result = (c1i, c1j, c2i, c2j)
        return result
