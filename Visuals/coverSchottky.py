from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import diskModel as dm
import poincareModel as pm
from mpmath import mp

'''
This file contains functions for visualizing the triple cover of a symmetric
Schottky group generated by 3 reflections. The cover is obtained by including
rotation by 2pi/3 as a group element.
'''

# computes inverse Cayley transform
def inv_Cayley(w) :
    if w == 1 :
        return np.inf
    return 1j*(1 + w)/(1 - w)

# maps fundamental domain in disk model to upper half plane
def disk_to_UHP(fd) :
    gs = []

    # apply inverse Cayley transform to each geodesic
    for g in fd.geodesics :
        gs.append(pm.Geodesic(inv_Cayley(g.ep1).real, inv_Cayley(g.ep2).real))

    # if fd has a reference point, map that to upper half plane
    ref_point = None
    if fd.referencePoint != None :
        ref_point = inv_Cayley(fd.referencePoint)

    return pm.FundamentalDomain(gs, ref_point)

# given a hyperbolic matrix A, returns the fixed points of the Mobius
# transformation defined by A
def get_fixed_points(A) :
    tr = A[0, 0] + A[1, 1]

    # matrix must be hyperbolic
    assert abs(tr) > 2, 'matrix must be hyperbolic'

    # if A[1, 0] is 0, infinity and -A[0, 1]/A[0, 0] are fixed points
    if A[1, 0] == 0 :
        return -A[0, 1]/A[0, 0], np.inf

    # otherwise, use quadratic formula to compute fixed points
    z1 = (A[0, 0] - A[1, 1] - np.sqrt(tr**2 - 4))/(2*A[1, 0])
    z2 = (A[0, 0] - A[1, 1] + np.sqrt(tr**2 - 4))/(2*A[1, 0])
    return z1, z2

# returns fundamental domain in the disk for the group generated by
#   - one circle from the symmetric Schottky group with 3 circles of angle
#     thet (we use the rightmost circle, symmetric across the real axis)
#   - two reflections across diameters which compose to obtain rotation by
#     2pi/3 (we use the real line and the rotation of that by pi/3)
def get_reflection_group(thet) :
    # form the geodesics
    R = dm.Geodesic(np.exp(1j*thet/2), np.exp(-1j*thet/2))
    D1 = dm.Geodesic(np.exp(1j*np.pi/3), -np.exp(1j*np.pi/3))
    D2 = dm.Geodesic(-1, 1)

    # return the fundamental domain with an appropriate reference point
    return dm.FundamentalDomain([R, D1, D2], 0.1 + 1j*0.1)

# draws a fundamental domain in the disk for the group as given by
# get_reflection_group()
def draw_reflection_group(thet) :
    # get fundamental domain
    fd = get_reflection_group(thet)

    # plot it!
    ax = dm.setupFig()
    fd.draw(ax)

    # label geodesics
    plt.text(np.cos(thet/2) + 0.02, np.sin(thet/2) + 0.02, '$R$', fontsize=12)
    plt.text(0.52, 0.88, '$D_1$', fontsize=12)
    plt.text(1.02, 0, '$D_2$', fontsize=12)

    plt.axis('off')
    plt.show()

# returns fundamental domain for group obtained by doubling reflection group
# across D1
def get_doubled_group(thet) :
    # form the geodesics for the fundamental domain
    R = dm.Geodesic(np.exp(1j*thet/2), np.exp(-1j*thet/2))
    D1R = dm.Geodesic(np.exp(1j*(2*np.pi/3 - thet/2)), np.exp(1j*(2*np.pi/3 + thet/2)))
    D2 = dm.Geodesic(-1, 1)
    D1D2 = dm.Geodesic(np.exp(1j*2*np.pi/3), -np.exp(1j*2*np.pi/3))

    # return the fundamental domain with an appropriate reference point
    return dm.FundamentalDomain([R, D1R, D2, D1D2], 0.6 + 1j*0.6)

# draws a fundamental domain in the disk for the group as given by
# get_doubled_group()
def draw_doubled_group(thet, show_axis=False) :
    # get fundamental domain
    fd = get_doubled_group(thet)

    # we'll also plot D1 as a dashed line
    D1 = dm.Geodesic(np.exp(1j*np.pi/3), -np.exp(1j*np.pi/3))

    # plot stuff!
    ax = dm.setupFig()
    fd.draw(ax)
    D1.draw(ax, lstyle='--')

    # if asked, show axis of D1R (this cuts off the flare)
    if show_axis :
        M = 0.5*np.array([ [(-np.sqrt(3) + 1j)/np.sin(thet/2), (np.sqrt(3) - 1j)/np.tan(thet/2)],
                           [(np.sqrt(3) + 1j)/np.tan(thet/2), -(np.sqrt(3) + 1j)/np.sin(thet/2)] ])
        z1, z2 = get_fixed_points(M)
        D = dm.Geodesic(z1, z2)
        D.draw(ax, lstyle=':')
        plt.text(0.98, 0.28, 'Axis of $D_1R$', fontsize=12)

    # label geodesics
    plt.text(np.cos(thet/2) + 0.02, np.sin(thet/2) + 0.02, '$R$', fontsize=12)
    plt.text(0.52, 0.88, '$D_1$', fontsize=12)
    plt.text(1.02, 0, '$D_2$', fontsize=12)
    plt.text(np.cos(2*np.pi/3 - thet/2), np.sin(2*np.pi/3 - thet/2) + 0.02, '$D_1(R)$', fontsize=12)
    plt.text(-0.56, 0.94, '$D_1(D_2)$', fontsize=12)

    plt.axis('off')
    plt.show()

# draws a fundamental domain in the disk for the group as given by
# get_doubled_group() after mapping to upper half plane
def doubled_group_UHP(thet, show_axis=False, X=4, Y=4) :
    # get fundamental domain
    fd = disk_to_UHP(get_doubled_group(thet))

    # we'll also plot image of D1 as a dashed line
    D1 = pm.Geodesic(inv_Cayley(np.exp(1j*np.pi/3)).real,
                     inv_Cayley(-np.exp(1j*np.pi/3)).real)

    # plot stuff!
    ax = pm.setupFig(-X, X, 0, Y)
    fd.draw(ax, fill=1)
    D1.draw(ax, lstyle='--')

    # label geodesics
    label_pts = []
    labels = []
    label_pts.append(inv_Cayley(np.cos(thet/2) + 1j*np.sin(thet/2)).real)
    labels.append('$R$')
    label_pts.append(inv_Cayley(0.5 + 1j*np.sqrt(3)/2).real)
    labels.append('$D_1$')
    label_pts.append(0)
    labels.append('$D_2$')
    label_pts.append(inv_Cayley(np.cos(2*np.pi/3 - thet/2) + 1j*np.sin(2*np.pi/3 - thet/2)).real)
    labels.append('$D_1(R)$')
    label_pts.append(inv_Cayley(0.5 - 1j*np.sqrt(3)/2).real)
    labels.append('$D_1(D_2)$')

    # if asked, show axis of D1R (this cuts off the flare)
    if show_axis :
        M = 0.5*np.array([ [-np.sqrt(3)*np.tan(thet/4), 1/np.tan(thet/4)],
                           [-np.tan(thet/4), -np.sqrt(3)/np.tan(thet/4)] ])
        z1, z2 = get_fixed_points(M)
        D = pm.Geodesic(z1, z2)
        D.draw(ax, lstyle=':')
        label_pts.append(z2)
        labels.append('Axis of $D_1R$')

    plt.xticks(label_pts, labels)
    plt.yticks([])

    return ax

# draws the flare domain for doubled group
def draw_flare_domain(thet, X, Y) :
    # get fundamental domain
    fd = disk_to_UHP(get_doubled_group(thet))

    # get endpoints of D1R
    M = 0.5*np.array([ [-np.sqrt(3)*np.tan(thet/4), 1/np.tan(thet/4)],
                       [-np.tan(thet/4), -np.sqrt(3)/np.tan(thet/4)] ])
    z2, z1 = get_fixed_points(M) # swapped names to match Tex document
    D = pm.Geodesic(z1, z2)

    # t is the negative endpoint of R
    t = (fd.geodesics[0]).ep1

    # Define Mobius transformation and apply to whole fundamental domain
    mult = (t - z2)/(t - z1)
    U = np.array([ [mult, -z1*mult], [1, -z2] ])
    fd_flare = pm.mobius(U, fd)

    # Draw the fundamental domain!
    ax = pm.setupFig(-X, X, 0, Y)
    fd_flare.draw(ax, fill=1)

    # Let's draw D1 and the axis as well
    D1 = pm.Geodesic(inv_Cayley(np.exp(1j*np.pi/3)).real,
                     inv_Cayley(-np.exp(1j*np.pi/3)).real)
    D1_flare = pm.mobius(U, D1)
    D1_flare.draw(ax, lstyle='--')
    D_flare = pm.mobius(U, D)
    D_flare.draw(ax, lstyle=':')

    # label everything (R, D1R, D2, D1D2)
    labels = ['$U(R)$', '$U(D_1R)$', '$U(D_2)$', '$U(D_1D_2)$', '$U(D_1)$']
    label_pts = []
    for i in range(4) :
        g = fd_flare.geodesics[i]
        if i == 1 or i == 2 :
            label_pts.append(min(g.ep1, g.ep2))
        else :
            label_pts.append(max(g.ep1, g.ep2))
    label_pts.append(max(D1_flare.ep1, D1_flare.ep2))

    plt.xticks(label_pts, labels)
    plt.yticks([])

    plt.show()

# plots points from testPoints.txt along with fundamental domain in the UHP
def plot_test_points(thet) :
    # read file with points
    with open('testPoints.txt') as f :
        line = f.readline().strip()

    # convert point strings to complex numbers
    pts_string = line.replace('I', '1j')
    pts_string = pts_string.replace('*', '').split(',')
    pts = []
    for p in pts_string :
        p_new = ''.join(p.split())
        pts.append(complex(p_new))

    # get pullbacks
    R, D1, _ = get_circles_UHP(thet)
    pts_pb = []
    for pt in pts :
        pts_pb.append(pullback_Schottky(pt, R, D1))

    # draw fundamental domain
    ax = doubled_group_UHP(thet)

    # draw points
    pts = np.array(pts)
    pts_pb = np.array(pts_pb)
    ax.plot(pts.real, pts.imag, '.k')
    ax.plot(pts_pb.real, pts_pb.imag, 'xk')

    plt.show()

###################################
# Functions for testing PARI code #
###################################

# computes reflection through geodesic semicircle in upper half plane
def reflect(z, cent, rad) :
    if abs(z - cent) < 1e-10 :
        return np.inf
    return rad**2 / (z.real - 1j*z.imag - cent) + cent

# get center and radius for R, D1, and D2 in UHP model
def get_circles_UHP(thet) :
    # get fundamental domain
    fd = disk_to_UHP(get_reflection_group(thet))

    # recall that geodesic order is R, D1, D2

    # get R
    g = fd.geodesics[0]
    R = [0, 0]
    R[0] = (g.ep1 + g.ep2)/2
    R[1] = abs(g.ep1 - g.ep2)/2

    # get D1
    g = fd.geodesics[1]
    D1 = [0, 0]
    D1[0] = (g.ep1 + g.ep2)/2
    D1[1] = abs(g.ep1 - g.ep2)/2

    # get D2
    D2 = [0, np.inf]

    return R, D1, D2

# try to find pullbacks of test points
def specific_pullbacks() :
    R, D1, _ = get_circles_UHP(np.pi/6)

    # interested in the following three points
    z1 = 1j/10
    z2 = 1/2 + 1j
    z3 = 1 + 1j

    # pullback for z1 is D2D1
    print 'D2D1(z1) ='
    D1z1 = reflect(z1, D1[0], D1[1])
    print str(-D1z1.real + 1j*D1z1.imag)
    print ''

    # pullback for z2 is D1D2
    print 'D1D2(z2)'
    print reflect(-z2.real + 1j*z2.imag, D1[0], D1[1])
    print ''

    # pullback for z3 is D1D2
    print 'D1D2(z3)'
    print reflect(-z3.real + 1j*z3.imag, D1[0], D1[1])
    print ''

# general pullback algorithm, which goes by
#   - if z is in D1, reflect by D1
#   - elif z is outside R, reflect by R
#   - elif Re(z) > 0, reflect across imaginary axis
#   - else, in original FD, so...
#       - if used odd number of reflections, return D1(z)
#       - else, return z
def pullback_Schottky(z, R, D1) :
    # save z at start in case we enter an infinite loop
    z_start = z

    # how many reflections have we performed so far?
    num_refl = 0

    # loop until done
    while 1 :

        # check if in D1
        if abs(z - D1[0]) < D1[1] :
            z = reflect(z, D1[0], D1[1])
            num_refl += 1

        # check if outside R
        elif abs(z - R[0]) > R[1] :
            z = reflect(z, R[0], R[1])
            num_refl += 1

        # check if real part positive
        elif z.real > 0 :
            z = -z.real + 1j*z.imag
            num_refl += 1

        # otherwise we are in original fundamental domain
        else :
            # if odd number of reflections, need to reflect by D1
            if num_refl % 2 == 1 :
                return reflect(z, D1[0], D1[1])
            else :
                return z

        # if we ever end up back at start, in an infinite loop
        if abs(z - z_start) < 1e-10 :
            print 'Warning, pullback algorithm for ' + str(z_start) + ' returned to self.'
            return z_start

# Given two (center, radius) pairs (c, r) and (a, t), computes the fixed
# points of the Mobius transformation obtained from reflection across the
# first circle composed with reflection across the second circle
def fixed_points_refls(c, r, a, t) :
    A = (c**2 - a**2 + t**2 - r**2)/2/(c - a)
    B = np.sqrt(( (a - c)**2 - (r**2 + t**2) )**2 - 4*r**2*t**2 )/2/(c - a)
    return A - B, A + B

# get flare parameters
#   - w1, w2: left- and right-hand points of axis cutting off the flare
#   - t: leftmost point of flare
#   - pre_kappa: rightmost point of flare
def get_flare_data(thet) :
    R, D1, _ = get_circles_UHP(np.pi/6)

    # axis is D1R
    z1, z2 = fixed_points_refls(D1[0], D1[1], R[0], R[1])
    w1 = min(z1, z2)
    w2 = max(z1, z2)

    # t is leftmost point of R
    t = R[0] - R[1]

    # pre_kappa is D1(t)
    pre_kappa = reflect(t, D1[0], D1[1])

    return w1, w2, t, pre_kappa

# map points to the flare expansion
#   - w1, w2: left- and right-hand points of axis cutting off the flare
#   - t: leftmost point of flare
def map_to_flare(z, w1, w2, t) :
    w = (t - w2)/(t - w1)*(z - w1)/(z - w2)
    return [abs(w), np.arccos(w.real/abs(w))]

# use mpmath to compute Whittaker function for the flare expansion
def flare_Whitt(thet, m, s, kappa) :
    nu = s - 1/2
    mu = -1/2 + 2*np.pi*1j*m/np.log(kappa)
    return np.sqrt(np.sin(thet))*mp.legenp(mu, -nu, np.cos(thet)).real

# approximates values which appear in test_init_eqns in
# Algorithms/test_cover_Schottky.pari
def test_init_eqns() :
    # points of interest and their (approximate) pullbacks
    z1 = 1j/10
    z2 = 1/2 + 1j
    z3 = 1 + 1j
    zpb1 = -1.6647866985370762335846328913503045080 + 0.38834951456310679611650485436893203883*1j
    zpb2 = -0.47482996250770193378102435128659914646 + 1.3254033600140835360567302423526384393*1j
    zpb3 = -1.0554745644123689251027543177449890821 + 1.1312542286635410364383667765418102415*1j

    # check these points in the flares
    w1, w2, t, pre_kappa = get_flare_data(np.pi/6)
    z1_flare = map_to_flare(z1, w1, w2, t)
    z2_flare = map_to_flare(z2, w1, w2, t)
    z3_flare = map_to_flare(z3, w1, w2, t)
    zpb1_flare = map_to_flare(zpb1, w1, w2, t)
    zpb2_flare = map_to_flare(zpb2, w1, w2, t)
    zpb3_flare = map_to_flare(zpb3, w1, w2, t)
    zs = [z1_flare, z2_flare, z3_flare]
    zpbs = [zpb1_flare, zpb2_flare, zpb3_flare]

    # print info
    print 'z1 and zpb1 in flare'
    print z1_flare
    print zpb1_flare
    print ''
    print 'z2 and zpb2 in flare'
    print z2_flare
    print zpb2_flare
    print ''
    print 'z3 and zpb3 in flare'
    print z3_flare
    print zpb3_flare
    print ''

    # now form the matrix for the linear system
    kappa = map_to_flare(pre_kappa, w1, w2, t)[0]
    A = np.zeros((3, 3))
    for n in range(3) :
        z = zs[n]
        zpb = zpbs[n]
        for m in range(3) :
            A[n, m] = (flare_Whitt(z[1], m, 0.4, kappa)*
                       np.cos(2*np.pi*m*np.log(z[0])/np.log(kappa)) -
                       flare_Whitt(zpb[1], m, 0.4, kappa)*
                       np.cos(2*np.pi*m*np.log(zpb[0])/np.log(kappa)))

    # print matrix
    print 'A ='
    for n in range(3) :
        print A[n, :]

if __name__ == '__main__' :
    plot_test_points(116*np.pi/180)
