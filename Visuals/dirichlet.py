import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
import poincareModel as pm

'''
This file contains functions for constructing a Dirichlet domain for groups
acting discretely on the upper half-plane model of hyperbolic 2-space. We
consider the explicit example given in
    https://sites.math.rutgers.edu/~alexk/files/UniformLattice.pdf
'''

#

# finds the symbolic expression for the endpoints of the geodesic which is
# equidistant from z1 and z2
def symbolic_equidistant(z1, z2) :
    # get names for real and imagniary parts for ease of reading
    x1, y1 = z1.real, z1.imag
    x2, y2 = z2.real, z2.imag

    # split into cases depending on if geodesic G through z1 and z2 is a line
    # or a half-circle

    # G is a line
    if x1 == x2 :

        # get point z0 on G equidistant from z1 and z2
        x0 = x1
        y0 = sp.sqrt(y1*y2)

        # z0 is at top of circle, so center is x0 and radius is y0
        return x0 - y0, x0 + y0

    # G is a half-circle
    else :

        # get point z0 on G equidistant from z1 and z2
