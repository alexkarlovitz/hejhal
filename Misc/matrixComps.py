from __future__ import division
import numpy as np
from sympy import *

'''
Some matrix computations which come up in this project.
'''

# make fraction in sympy
def R(a, b) :
    return Rational(a, b)

# matrix for the Cayley transform
def get_C() :
    return Matrix( [[1, -I],
                    [1, I]] )

# matrix for the inverse Cayley transform
def get_Cinv() :
    return Matrix( [[R(1, 2), R(1, 2)],
                    [I*R(1, 2), -I*R(1, 2)]] )

# computes Mobius transformation A(z)
def mobiusTransform(A, z) :
    if z == np.inf :
        if A[1, 0] == 0.0 :
            return np.inf
        return A[0, 0] / A[1, 0]
    elif A[1, 0] * z + A[1, 1] == 0.0 :
        return np.inf
    return (A[0, 0] * z + A[0, 1]) / (A[1, 0] * z + A[1, 1])

########################
# Schottky group cover #
########################

# my formula for D1D2
def get_D1D2() :
    return Matrix( [[R(1, 2) + I*sqrt(3)/2, 0],
                    [0, R(1, 2) - I*sqrt(3)/2]] )

# my formula for D2R
def get_D2R() :
    T = symbols('T')
    return Matrix( [[I*csc(T/2), -I*cot(T/2)],
                    [I*cot(T/2), -I*csc(T/2)]] )

# check conjugates to SL(2, R)
def generators_H() :
    D1D2 = get_D1D2()
    D2R = get_D2R()
    C = get_C()
    Cinv = get_Cinv()

    print 'D1D2'
    pprint(D1D2)
    print ''

    print 'CinvD1D2C'
    pprint(simplify(Cinv*D1D2*C))
    print ''

    print 'D2R'
    pprint(D2R)
    print ''

    print 'CinvD2RC'
    pprint(simplify(Cinv*D2R*C))

# what does D1R look like? (This is the hyperbolic matrix whose axis
# cuts off the flare)
def compute_D1R() :
    # since D2 is an involution, D1R = D1D2D2R
    D1D2 = get_D1D2()
    D2R = get_D2R()
    D1R = D1D2*D2R

    return simplify(D1R)

def hyperbolic_info() :
    D1R = compute_D1R()
    C = get_C()
    Cinv = get_Cinv()

    print 'D1R'
    pprint(D1R)
    print ''

    print 'In H:'
    pprint(simplify(Cinv*D1R*C))

if __name__ == '__main__' :
    hyperbolic_info()
