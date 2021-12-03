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

# applies the Mobius transformation A to the point z and gives symbolic solution
def symbolic_Mobius(A, z) :
    if z == np.inf :
        if A[1, 0] == 0.0 :
            return np.inf
        return A[0, 0] / A[1, 0]
    elif A[1, 0] * z + A[1, 1] == 0 :
        return np.inf
    return (A[0, 0] * z + A[0, 1]) / (A[1, 0] * z + A[1, 1])

# finds the symbolic expression for the endpoints of the geodesic which is
# equidistant from z1 and z2
def symbolic_equidistant(z1, z2) :
    # get names for real and imagniary parts for ease of reading
    x1, y1 = sp.re(z1), sp.im(z1)
    x2, y2 = sp.re(z2), sp.im(z2)

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

        # get center and radius of G
        c = (x2**2 + y2**2 - x1**2 - y1**2)/2/(x2 - x1)
        r = sp.sqrt( (x1 - c)**2 + y1**2 )

        # get point z0 on G equidistant from z1 and z2
        x = sp.Symbol('x')
        y = sp.Symbol('y')
        expr1 = y2*((x - x1)**2 + (y - y1)**2) - y1*((x - x2)**2 + (y - y2)**2)
        expr2 = (x - c)**2 + y**2 - r**2

        # the two circles will intersect at two points; we want the one with
        # a positive y
        solns = sp.solve([expr1, expr2])
        x0 = solns[0][x]
        y0 = max(solns[0][y], solns[1][y])

        # next, get the slope of the desired geodesic at z0
        s = (c - x0)/y0

        # if slope is 0, solution is a line (second endpoint infty)
        if s == 0:
            return x0, np.inf

        # otherwise, solve for center and radius
        a = -y0/s + x0
        R = sp.sqrt( (x0 - a)**2 + y0**2 )

        return a - R, a + R

#####################
# Cocompact Example #
#####################

# convert four real numbers into a sympy matrix
#   [ [alpha, -beta],
#     [bar(beta), bar(alpha)] ]
# where alpha = a + b*sqrt(3), beta = d + c*sqrt(3)
def vec_to_M(a, b, c, d) :
    return sp.Matrix([ [a + b*sp.sqrt(3), -d - c*sp.sqrt(3)],
                       [d - c*sp.sqrt(3), a - b*sp.sqrt(3)]])

# the list of matrices given in Kontorovich's example
#   - first matrix is S
#   - second is scaling matrix
def get_all_matrices() :
    L = [(0, 0, 0, 1),
         (2, 1, 0, 0),
         (2, -1, 0, 0),
         (-3, 0, -2, -2),
         (-3, 0, 2, 2),
         (-2, 0, -1, 0),
         (-2, 0, 1, 0),
         (-2, 0, -2, -3),
         (-2, 0, 2, 3),
         (-3, -2, 0, 2),
         (-3, -2, 0, -2),
         (0, 0, 1, 2)]

    Ms = []
    for v in L :
        Ms.append(vec_to_M(*v))

    return Ms

# get all geodesics equidistant from 2i and each point in its orbit under the
# matrices
def get_all_geodesics() :
    # get matrices
    Ms = get_all_matrices()

    # for each M, get geodesic eq from 2i and M(2i)
    gs = []
    z = 2*sp.I
    for M in Ms :
        Mz = symbolic_Mobius(M, z)
        gs.append(symbolic_equidistant(z, Mz))

    return gs

# draw every geodesic obtained from get_all_geodesics
def draw_all_geodesics() :
    # get all geodesics
    gs = get_all_geodesics()

    # set up a figure to plot on
    ax = pm.setupFig(-8, 8, 0, 8)

    # for each pair of points, create Geodesic object and draw on plot
    ticks = []
    labels = []
    for i in range(len(gs)) :
        g = gs[i]
        G = pm.Geodesic(float(sp.N(g[0])), float(sp.N(g[1])))
        G.draw(ax)

        # label x-axis with which matrix this is
        ticks.append(float(sp.N(g[0])))
        labels.append('$g_{' + str(i) + '}$')

    plt.xticks(ticks, labels)

    plt.show()

# get only the matrices which contribute to fundamental domain in Kontorovich's
# example
def get_fd_matrices() :
    L = [(0, 0, 0, 1),
         (-3, 0, -2, -2),
         (-3, 0, 2, 2),
         (-2, 0, -1, 0),
         (-2, 0, 1, 0),
         (-2, 0, -2, -3),
         (-2, 0, 2, 3),
         (0, 0, 1, 2)]

    Ms = []
    for v in L :
        Ms.append(vec_to_M(*v))

    return Ms

# get all geodesics equidistant from 2i and each point in its orbit under the
# matrices from get_fd_matrices
def get_fd_geodesics() :
    # get matrices
    Ms = get_fd_matrices()

    # for each M, get geodesic eq from 2i and M(2i)
    gs = []
    z = 2*sp.I
    for M in Ms :
        Mz = symbolic_Mobius(M, z)
        gs.append(symbolic_equidistant(z, Mz))

    return gs

# gets list of centers and radii of the geodesics from get_fd_geodesics
def get_fd_data() :
    # get geodesic end points
    gs = get_fd_geodesics()

    # get centers and radii
    cs = []
    rs = []
    for g in gs :
        cs.append(sp.simplify( (g[0] + g[1])/2 ))
        rs.append(sp.simplify( abs(g[1] - g[0])/2 ))

    return cs, rs

# draw only the matrices from get_fd_matrices
def draw_fd_geodesics() :
    # get all geodesics
    gs = get_fd_geodesics()

    # set up a figure to plot on
    ax = pm.setupFig(-6, 6, 0, 5)

    # for each pair of points, create Geodesic object and draw on plot
    ticks = []
    labels = []
    for i in range(len(gs)) :
        g = gs[i]
        G = pm.Geodesic(float(sp.N(g[0])), float(sp.N(g[1])))
        G.draw(ax)

        # label x-axis with which matrix this is
        ticks.append(float(sp.N(g[0])))
        labels.append('$g_{' + str(i) + '}$')

    plt.xticks(ticks, labels)

    plt.show()

# create a FundamentalDomain object for this example and draw it
def draw_fd() :

    # get all geodesics as Geodesic objects and keep track of x ticks
    pairs = get_fd_geodesics()
    gs = []
    ticks = []
    labels = []
    for i in range(len(pairs)) :
        p = pairs[i]
        gs.append(pm.Geodesic(float(sp.N(p[0])), float(sp.N(p[1]))))

        # label x-axis with which matrix this is
        ticks.append(float(sp.N(p[0])))
        labels.append('$g_{' + str(i) + '}$')

    # place in FundamentalDomain object with 2i as reference point
    fd = pm.FundamentalDomain(gs, 2j)

    # set up a figure to plot on
    ax = pm.setupFig(-6, 6, 0, 5)

    # draw and show the image
    fd.draw(ax, fill=1)
    plt.xticks(ticks, labels)

    return ax

# given an ordering of the Ms (with corresponding centers and radii of the
# geodesics), computes pullback of z to fundamental domain
def pullback(z, ordering, inverses, Ms, cs, rs, idx_7=1) :
    # keep track of actions
    word = ''

    # for clarity, reorder everything first
    cs_base = []
    rs_base = []
    for idx in ordering :
        cs_base.append(cs[idx])
        rs_base.append(rs[idx])

    Ms_inv = []
    for idx in inverses :
        Ms_inv.append(Ms[idx])

    # loop until done
    done = False
    while not done :

        # loop through ordering
        done = True
        for i in range(len(Ms)) :
            # check if in this geodesic
            in_geodesic = False
            if i == idx_7 :
                in_geodesic = abs(z - cs_base[i]) > rs_base[i]
            else :
                in_geodesic = abs(z - cs_base[i]) < rs_base[i]

            # if in this geodesic, map by inverse, then restart for loop
            if in_geodesic :
                done = False
                z = pm.mobiusTransform(Ms_inv[i], z)
                word = str(inverses[i]) + word
                break

    # if we got through the for loop without breaking, then done
    return z, word

#########
# Tests #
#########

def test_eq() :
    print(symbolic_equidistant(sp.I, sp.exp(2)*sp.I))
    print(symbolic_equidistant(sp.I, 1 + sp.I))
    print(symbolic_equidistant(sp.I, 2 + sp.Rational(1, 2)*sp.I))

def test_data() :
    Ms = get_all_matrices()
    for M in Ms :
        print(M)
    print()

    gs = get_all_geodesics()
    for g in gs :
        print(sp.N(g[0]), sp.N(g[1]))
    print()

    Ms = get_fd_matrices()
    for M in Ms :
        print(M)
    print()

    gs = get_fd_geodesics()
    for g in gs :
        print(sp.N(g[0]), sp.N(g[1]))
    print()

    cs, rs = get_fd_data()
    for i in range(len(cs)) :
        print(sp.N(cs[i]), sp.N(rs[i]))

def test_pullback() :
    # set up ordering, Ms, centers, and radii
    ordering = [0, 7, 3, 1, 5, 4, 2, 6]
    inverses = [0, 7, 4, 2, 6, 3, 1, 5]
    Ms = np.array([sp.N(M) for M in get_fd_matrices()], dtype=np.float)
    cs, rs = get_fd_data()
    cs = [float(sp.N(c)) for c in cs]
    rs = [float(sp.N(r)) for r in rs]

    # test z values
    kappa = (2 + np.sqrt(3))/(2 - np.sqrt(3))
    zs = [1j/2, 5.5 + 1j, kappa**(1/50)*np.exp(1j*np.pi/20)]

    # for each z, print pullback and word to get there
    for z in zs :
        zpb, word = pullback(z, ordering, inverses, Ms, cs, rs)
        print(z, zpb, word)
        print()

def draw_pullbacks() :
    # set up ordering, Ms, centers, and radii
    ordering = [0, 7, 3, 1, 5, 4, 2, 6]
    inverses = [0, 7, 4, 2, 6, 3, 1, 5]
    Ms = np.array([sp.N(M) for M in get_fd_matrices()], dtype=np.float)
    cs, rs = get_fd_data()
    cs = [float(sp.N(c)) for c in cs]
    rs = [float(sp.N(r)) for r in rs]

    # get N z values along ray of angle alpha evenly spaced in log space from
    # length 1 to sqrt(kappa)
    kappa = (2 + np.sqrt(3))/(2 - np.sqrt(3))

    thet = np.pi/20
    N = 50

    rhos = kappa**(np.arange(N)/(2*N))
    zs = [rho*np.exp(1j*thet) for rho in rhos]

    # get pullbacks
    zpbs = []
    for z in zs :
        zpb, word = pullback(z, ordering, inverses, Ms, cs, rs)
        print(word)
        zpbs.append(zpb)

    # plot fundamental domain and the points
    ax = draw_fd()
    zs = np.array(zs)
    zpbs = np.array(zpbs)

    ax.plot(zs.real, zs.imag, '.')
    ax.plot(zpbs.real, zpbs.imag, 'x')

    plt.show()
