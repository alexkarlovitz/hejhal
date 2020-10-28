from __future__ import division
import numpy as np
import sympy as sp
from mpmath import mp

########################################################################################
# File contains functions to verify formulas in diskApproach.pdf and flareApproach.pdf #
########################################################################################

##############################
# Funcs for disk Formula (1) #
##############################

# converts rectangular (x, y) coordinates to polar (r, theta) coordinates
def rectToPol(x, y) :
    r = np.sqrt(x**2 + y**2)
    theta = np.arccos(x/r)
    if y < 0 :
        theta = 2*np.pi - theta
    return r, theta

# converts polar (r, theta) coordinates to rectangular (x, y) coordinates
def polToRect(r, theta) :
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    return x, y

# Cayley transform to go from upper half plane model to disk model
def cayley(x, y) :
    z = x + 1j*y
    tau = (z - 1j) / (z + 1j)
    return tau.real, tau.imag

# Inverse of the Cayley transform
def invCayley(u, v) :
    tau = u + 1j*v
    z = 1j*(1 + tau) / (1 - tau)
    return z.real, z.imag

# Map from upper half plane model to polar coordinates in disk model
def upperToPolar(x, y) :
    u, v = cayley(x, y)
    return rectToPol(u, v)

# Map from polar coordinates in disk model to upper half plane model
def polarToUpper(r, theta) :
    u, v = polToRect(r, theta)
    return invCayley(u, v)

# Conjectured formula from (r, theta) coordinates to (x, y) coordinates
def conjPTU(r, theta) :
    x = (-2*r*np.sin(theta)) / (1 + r**2 - 2*r*np.cos(theta))
    y = (1 - r**2) / (1 + r**2 - 2*r*np.cos(theta))
    return x, y

# Conjectured formula from (x, y) coordinates to (r, theta) coordinates
def conjUTP(x, y) :
    r = np.sqrt( (x**2 + (y-1)**2) / (x**2 + (y+1)**2) )
    if x**2 + y**2 == 1 :
        theta = np.pi / 2
        if x < 0 :
            theta = -theta
        return r, theta

    tanTheta = (-2*x) / (x**2 + y**2 - 1)
    theta = np.arctan(tanTheta)
    if (x**2 + y**2 - 1) < 0 :
        theta = theta + np.pi
    theta = theta % (2*np.pi)
    return r, theta

##############################
# Funcs for disk Formula (2) #
##############################
def recursionOutsideSum(r, n, L, c0, c1, c2, c3) :
    x1 = 2*c2*r**2 + 6*c2*r**3
    x2 = c1*r + 2*c2*r**2 + 3*c3*r**3
    x3 = -2*c1*r**3
    x4 = -n**2*(c0 + c1*r + c2*r**2 + c3*r**3)
    x5 = 2*n**2*(c0*r**2 + c1*r**3)
    x6 = -4*L*(c0*r**2 + c1*r**3)
    return x1 + x2 + x3 + x4 + x5 + x6

def recursionInsideSum(n, m, L, cm, cm_2, cm_4) :
    y1 = m*(m-1)*cm - 2*(m-2)*(m-3)*cm_2 + (m-4)*(m-5)*cm_4
    y2 = m*cm - 2*(m-2)*cm_2 + (m-4)*cm_4
    y3 = n**2*(-cm + 2*cm_2 - cm_4) - 4*L*cm_2
    return y1 + y2 + y3

def pw(a, b) :
    return mp.power(a, b)

def p(r, n, s) :
    return pw(1 - pw(r, 2), s) * pw(r, n)

def p_r(r, n, s) :
    return -2*s*pw(1-pw(r,2), s-1)*pw(r, n+1) + n*pw(1 - pw(r, 2), s) * pw(r, n-1)

def p_rr(r, n, s) :
    return 4*s*(s-1)*pw(1-pw(r, 2), s-2)*pw(r, n+2) - 2*(2*n+1)*s*pw(1-pw(r, 2), s-1)*pw(r, n) + n*(n-1)*pw(1-pw(r, 2), s) * pw(r, n-2)

def F(r, n, s) :
    return mp.hyp2f1(s, s+n, 1+n, pw(r, 2))

def F_r(r, n, s) :
    return 2*r*s*(s+n) / (1+n) * mp.hyp2f1(s+1, s+n+1, 2+n, pw(r, 2))

def F_rr(r, n, s) :
    return 2*s*(s+n) / (1+n) * mp.hyp2f1(s+1, s+n+1, 2+n, pw(r, 2)) + 4*pw(r, 2)*s*(s+1)*(s+n)*(s+n+1) / (1+n) / (2+n) * mp.hyp2f1(s+2, s+n+2, 3+n, pw(r, 2))

# Conjectured a_n numerical
def an(r, n, s) :
    return p(r, n, s) * F(r, n, s)

# Conjectured a_n' numerical
def an_r(r, n, s) :
    return p_r(r, n, s)*F(r, n, s) + p(r, n, s)*F_r(r, n, s)

# Conjectured a_n'' numerical
def an_rr(r, n, s) :
    return p_rr(r, n, s)*F(r, n, s) + 2*p_r(r, n, s)*F_r(r, n, s) + p(r, n, s)*F_rr(r, n, s)

# DE LHS
def deLHS(r, n, s) :
    return pw(1 - pw(r, 2), 2) / 4 *( an_rr(r, n, s) + an_r(r, n, s)/r - n**2*an(r, n, s)/pw(r, 2) )

# DE RHS
def deRHS(r, n, s) :
    return s*(s-1)*an(r, n, s)

# compare DE LHS and RHS
def compare_disk(r, n, s) :
    lhs = deLHS(r, n, s)
    rhs = deRHS(r, n, s)
    print lhs
    print rhs
    print abs(lhs - rhs)

###############################
# Funcs for flare Formula (2) #
###############################

# We write g_n = C*f(theta)*F(theta) where C is the part not depending on theta
# (doesn't matter in checking the diff eq), F is the hypergeometric function
# which appears, and f is the rest

# Use sympy to set up f and its derivatives
#   - using the symbol t for theta and M for mu/2
def get_f() :
    t, M = sp.symbols('t M')
    f = sp.sqrt(sp.sin(t)) * (1 + sp.cos(t))**M / (1 - sp.cos(t))**M
    f_t = sp.simplify( sp.diff(f, t) )
    f_tt = sp.simplify( sp.diff(f_t, t) )
    return f, f_t, f_tt

# Input tVal and mu into a sympy function of t and M = mu/2
def eval(func, t, tVal, M, mu) :
    return func.subs( [(t, tVal), (M, mu/2)] )

# Evaluate F at given t, nu, mu
def F2(t, nu, mu) :
    return mp.hyp2f1(-nu, nu+1, 1 - mu, (1 - mp.cos(t))/2)

# Evaluate F' at given t, nu, mu
def F2_t(t, nu, mu) :
    return -nu * (nu+1) / 2 / (1 - mu) * mp.hyp2f1(-nu+1, nu+2, 2-mu, (1 - mp.cos(t))/2) * mp.sin(t)

# Evaluate F'' at given t, nu, mu
def F2_tt(t, nu, mu) :
    return (-nu * (nu+1) / 2 / (1 - mu) * mp.hyp2f1(-nu+1, nu+2, 2-mu, (1 - mp.cos(t))/2) * mp.cos(t)
    - nu * (nu+1) * (-nu+1) * (nu+2) / 4 / (1-mu) / (2-mu) * mp.hyp2f1(-nu+2, nu+3,
    3-mu, (1 - mp.cos(t))/2) * mp.sin(t)**2 )

# DE (returns LHS and RHS)
def de(tVal, nu, mu, s, n, k) :
    # get f and its derivs from sympy
    f, f_t, f_tt = get_f()

    # evaluate f and its derivatives at our inputs
    t, M = sp.symbols('t M')
    fVal = eval(f, t, tVal, M, mu)
    f_tVal = eval(f_t, t, tVal, M, mu)
    f_ttVal = eval(f_tt, t, tVal, M, mu)

    # evaluate F and its derivatives at our inputs
    FVal = F2(tVal, nu, mu)
    F_tVal = F2_t(tVal, nu, mu)
    F_ttVal = F2_tt(tVal, nu, mu)

    # g(t) is f(t)*F(t)
    g = fVal * FVal
    g_tt = f_ttVal*FVal + 2*f_tVal*F_tVal + fVal*F_ttVal

    # evaluate conjectured expression
    lhs = mp.sin(tVal)**2*(g_tt + (2*mp.pi*1j*n / mp.log(k))**2*g)
    rhs = s*(s-1)*g
    return lhs, rhs

# main function
if __name__ == '__main__' :
    # set up test values
    r = 0.35
    s = 0.76

    # Test!
    for n in range(10) :
        print "n =", n
        compare_disk(r, n, s)
        print " "
