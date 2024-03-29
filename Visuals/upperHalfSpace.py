from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt

from sympy.abc import i, j, k

###########
# Classes #
###########

class Quaternion :
    """Quaternion algebra"""

    # the quaternion is x1 + x2i + yj
    def __init__(self, x1=0, x2=0, y=0, z=0) :
        self.x1 = x1
        self.x2 = x2
        self.y = y
        self.z = z

    # quaternions are equal if all components are
    def __eq__(self, other) :
        if not isinstance(self, Quaternion) or not isinstance(other, Quaternion) :
            return False
        if self.x1 != other.x1 or self.x2 != other.x2 or self.y != other.y or self.z != other.z :
            return False

        return True

    def __ne__(self, other) :
        return not self == other

    # get norm of quaternion
    def norm(self) :
        return np.sqrt(self.x1**2 + self.x2**2 + self.y**2 + self.z**2)

    # overload string
    def __str__(self) :
        # form sympy object
        q = self.x1 + self.x2*i + self.y*j + self.z*k

        return str(q)

    # overload addition
    def __add__(self, other) :
        return Quaternion(self.x1 + other.x1,
                          self.x2 + other.x2,
                          self.y + other.y,
                          self.z + other.z)

    # overload subtraction
    def __sub__(self, other) :
        return Quaternion(self.x1 - other.x1,
                          self.x2 - other.x2,
                          self.y - other.y,
                          self.z - other.z)

    # overload negation
    def __neg__(self) :
        return Quaternion(-self.x1, -self.x2,
                          -self.y, -self.z)

    # overload multiplication
    def __mul__(self, other) :
        return Quaternion(self.x1*other.x1 - self.x2*other.x2 - self.y*other.y - self.z*other.z,
                          self.x1*other.x2 + self.x2*other.x1 + self.y*other.z - self.z*other.y,
                          self.x1*other.y - self.x2*other.z + self.y*other.x1 + self.z*other.x2,
                          self.x1*other.z + self.x2*other.y - self.y*other.x2 + self.z*other.x1)

    # overload division
    def __truediv__(self, other) :
        D = other.norm()**2
        inv_other = Quaternion(other.x1/D, -other.x2/D, -other.y/D, -other.z/D)
        return self*inv_other

class GeodesicWall :
    """A single geodesic wall in the upper half space model of hyperbolic 3 space."""

    # initialize a geodesic wall by specifying three points on the x1x2-plane;
    # if points are on a line or one point is infinity, the geodesic wall is
    # a half plane; otherwise, it is a half sphere
    def __init__(self, q1, q2, q3) :

        # check that points are distinct
        if q1 == q2 or q1 == q3 or q2 == q3:
            raise ValueError('Points defining a geodesic wall must be distinct.')

        # if point isn't infinity, it must be a j- and k-free quaternion
        self.verify_point(q1)
        self.verify_point(q2)
        self.verify_point(q3)

        # initialize variables
        self.q1 = q1
        self.q2 = q2
        self.q3 = q3

    # checks if point is either inf or a j- and k-free quaternion
    def verify_point(self, q) :
        if q != np.inf :
            if not isinstance(q, Quaternion) :
                raise ValueError('Point ' + str(q) + ' is neither infinity nor a Quaternion.')

            if q.y !=0 or q.z != 0 :
                raise ValueError('Points defining a geodesic wall must be j- and k-free. Received point ' + str(q) + '.')

    # returns [center, radius, True] where [center, radius] is of circle
    # describing this hemisphere; if it's a vertical half plane, returns
    # [(x1, y1), (x2, y2), False] where (x1, y1) and (x2, y2) are
    # non-infinite points on the line
    #   see http://web.archive.org/web/20161011113446/http://www.abecedarical.com/zenosamples/zs_circle3pts.html
    def get_geometric_desc(self) :
        # check if one of the points is inf
        if self.q1 == np.inf or self.q2 == np.inf or self.q3 == np.inf :
            # find non-infinite points
            p1, p2 = self.q1, self.q2
            if p1 == np.inf :
                p1 = self.q3
            elif p2 == np.inf :
                p2 = self.q3

            # return non-infinite points
            return [(p1.x1, p1.x2), (p2.x1, p2.x2), False]

        # otherwise, all three are j- and k-free Quaternions
        else :
            # write points as pairs (x, y) for ease of reading
            x1, y1 = self.q1.x1, self.q1.x2
            x2, y2 = self.q2.x1, self.q2.x2
            x3, y3 = self.q3.x1, self.q3.x2

            # first cofactor determines line or circle
            M11 = np.linalg.det([ [x1, y1, 1],
                                  [x2, y2, 1],
                                  [x3, y3, 1] ])

            # if M11 = 0, the points are in a line
            if M11 == 0 :
                return [(x1, y1), (x2, y2), False]

            # otherwise, points are in a circle
            else :
                M12 = np.linalg.det([ [x1**2 + y1**2, y1, 1],
                                      [x2**2 + y2**2, y2, 1],
                                      [x3**2 + y3**2, y3, 1] ])
                M13 = np.linalg.det([ [x1**2 + y1**2, x1, 1],
                                      [x2**2 + y2**2, x2, 1],
                                      [x3**2 + y3**2, x3, 1] ])
                M14 = np.linalg.det([ [x1**2 + y1**2, x1, y1],
                                      [x2**2 + y2**2, x2, y2],
                                      [x3**2 + y3**2, x3, y3] ])

                # apply formulas for center and radius
                x = M12/M11/2
                y = -M13/M11/2
                r = np.sqrt(x**2 + y**2 + M14/M11)

                # return
                return [(x, y), r, True]

    # draw this geodesic wall on a given matplotlib axis
    def draw(self, ax, color='k') :
        # get geometric description
        g = self.get_geometric_desc()

        # third element of g tells if circle
        is_circle = g[2]

        # draw either half sphere or half plane
        if is_circle :
            self.draw_sphere(ax, g[0][0], g[0][1], g[1], color)
        else :
            self.draw_plane(ax, g[0][0], g[0][1], g[1][0], g[1][1], color)

    # given two points on x1x2-plane, draw vertical half plane through them
    def draw_plane(self, ax, x1, y1, x2, y2, color='k') :
        # get bounds for cube
        x_min, x_max = ax.get_xlim()
        y_min, y_max = ax.get_ylim()
        z_min, z_max = ax.get_zlim()

        # check if plane is parallel to x1 axis
        if y1 == y2 :
            # make mesh in x and z coordinates
            X, Z = np.mgrid[x_min:x_max:2j, z_min:z_max:2j]

            # now plot
            Y = np.zeros_like(X) + y1
            ax.plot_surface(X, Y, Z, color=color, alpha=0.2)

        # check if plane is parallel to x2 axis
        elif x1 == x2 :
            # make mesh in y and z coordinates
            Y, Z = np.mgrid[y_min:y_max:2j, z_min:z_max:2j]

            # now plot
            X = np.zeros_like(Y) + x1
            ax.plot_surface(X, Y, Z, color=color, alpha=0.2)

        # otherwise, treat y as a function of x
        else :
            # make mesh in x and z coordinates
            X, Z = np.mgrid[x_min:x_max:2j, z_min:z_max:2j]

            # compute y and plot
            m = (y2 - y1)/(x2 - x1)
            Y = m*(X - x1) + y1
            ax.plot_surface(X, Y, Z, color=color, alpha=0.2)

    # given center (on x1x2-plane) and radius, draw half sphere
    def draw_sphere(self, ax, x, y, r, color='k') :
        # form meshgrid for unit half sphere, then scale and shift
        u, v = np.mgrid[0:2*np.pi:51j, 0:np.pi/2:51j]
        X = r*np.cos(u)*np.sin(v) + x
        Y = r*np.sin(u)*np.sin(v) + y
        Z = r*np.cos(v)

        # plot it!
        ax.plot_surface(X, Y, Z, color=color, alpha=0.2)

    # return string with the endpoints of the geodesic
    def __str__(self) :
        return str(self.q1) + ', ' + str(self.q2) + ', ' + str(self.q3)

#####################
# Mobius Transforms #
#####################

# given invertible 2x2 matrix A, returns A(z) where action is by Mobius transformation
#   - A has entries in C
#   - z is a Quaternion or np.inf
def mobiusTransform(A, z, min_size=1e-14) :
    # convert entries in A to Quaternions
    B = [[Quaternion(A[0, 0].real, A[0, 0].imag), Quaternion(A[0, 1].real, A[0, 1].imag)],
         [Quaternion(A[1, 0].real, A[1, 0].imag), Quaternion(A[1, 1].real, A[1, 1].imag)]]

    if z == np.inf :
        if A[1, 0] == 0.0 :
            return np.inf
        return B[0][0] / B[1][0]
    elif (B[1][0] * z + B[1][1]).norm() < min_size :
        return np.inf
    return (B[0][0] * z + B[0][1]) / (B[1][0] * z + B[1][1])

# given invertible 2x2 matrix A, returns mobius transformation A of a geodesic
# wall or list of such
def mobius(A, x) :
    # if x is a geodesic wall, get new geodesic where we applied A to endpoints
    if isinstance(x, GeodesicWall) :
        q1 = mobiusTransform(A, x.q1)
        q2 = mobiusTransform(A, x.q2)
        q3 = mobiusTransform(A, x.q3)
        return GeodesicWall(q1, q2, q3)

    # if x is a list of walls, apply A to every wall in list
    elif hasattr(x, '__iter__') :
        l = []
        for w in x :
            l.append(mobius(A, w))
        return l

    else :
        raise ValueError('Function "mobius" expects second input to be a GeodesicWall object (or a list of such). Received: ' + str(x))

# given invertible 2x2 matrix A, applies mobius transformation A to a list of
# geodesics

###############
# Reflections #
###############

# computes reflection of the Quaternion z through a GeodesicWall
def reflect(z, W, min_size=1e-14) :
    if not isinstance(z, Quaternion) :
        raise ValueError('First argument ' + str(z) + ' to "reflect" is not a Quaternion.')
    if not isinstance(W, GeodesicWall) :
        raise ValueError('Second argument ' + str(W) + ' to "reflect" is not a GeodesicWall.')

    # get geometric info about W
    g = W.get_geometric_desc()

    # proceed depending on whether W is a half sphere or half plane
    is_sphere = g[2]
    if is_sphere :
        if (z - c).norm() < min_size :
            return np.inf

        c = Quaternion(g[0][0], g[0][1])
        r_squared = Quaternion(g[1]**2)

        return r_squared*(z - c)/Quaternion((z - c).norm()**2) + c

    else :
        # get coordinates of point and points on line for ease of reading
        x, y = z.x1, z.x2
        x1, y1 = g[0][0], g[0][1]
        x2, y2 = g[1][0], g[1][1]

        # get vector from a point on line to point and vector on line
        ax, ay = x - x1, y - y1
        Lx, Ly = x2 - x1, y2 - y1

        # get orthogonal projection then orthogonal component
        p_scale = (ax*Lx + ay*Ly)/(Lx**2 + Ly**2)
        px, py = Lx*p_scale, Ly*p_scale
        ox, oy = ax - px, ay - py

        # reflection of z is z minus twice the orthogonal component
        rx, ry = zx - 2*ox, zy - 2*oy
        return Quaternion(rx, ry, z.y)

# computes reflection of the GeodesicWall W1 through the GeodesicWall W2
def reflect_wall(W1, W2) :
    g = W2.get_geometric_desc()
    is_sphere = g[2]

    # reflect each base point in W1
    qs = [W1.q1, W2.q2, W3.q3]
    new_qs = []
    for q in qs :
        # if infinity, can get answer immediately
        if q == np.inf :
            if is_sphere :
                new_qs.append(Quaternion(g[0][0], g[0][1]))
            else :
                new_qs.append(np.inf)

        # otherwise, use reflect function
        else :
            new_qs.append(reflect(q, W2))

    return GeodesicWall(new_qs[0], new_qs[1], new_qs[2])

############
# Plotting #
############

# prepare axis for 3d images
def axis_3d(x_lims=(-1, 1), y_lims=(-1, 1), z_lims=(-1, 1)) :
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.set_zlabel('$y$')
    ax.set_xlim(x_lims[0], x_lims[1])
    ax.set_ylim(y_lims[0], y_lims[1])
    ax.set_zlim(z_lims[0], z_lims[1])
    return ax

# plot list of GeodesicWalls
def plot_walls(Ws, b=None, x_lims=(-1, 1), y_lims=(-1, 1), z_lims=(-1, 1)) :
    # prepare 3d axis
    ax = None
    if b is None :
        ax = axis_3d(x_lims, y_lims, z_lims)
    else :
        ax = axis_3d((-b, b), (-b, b), (-1, 2*b - 1))

    # draw every wall
    for W in Ws :
        W.draw(ax)

    plt.show()

# plot list of walls on a given axis
def plot_Ws(ax, Ws, color='k') :
    for W in Ws :
        W.draw(ax, color)

# plot list of Quaternions
def plot_points(ax, qs, marker='o', color='k') :
    # make lists of x1, x2, and y-components
    xs = np.zeros(len(qs))
    ys = np.zeros(len(qs))
    zs = np.zeros(len(qs))
    for i in range(len(qs)) :
        q = qs[i]
        xs[i] = q.x1
        ys[i] = q.x2
        zs[i] = q.y

    ax.scatter(xs, ys, zs, marker=marker, color=color)

############
# Examples #
############

def example_fundamental_domain() :
    x1 = Quaternion(1, 0, 0, 0)
    x2 = Quaternion(0, 1, 0, 0)

    R = Quaternion(np.sqrt(2), 0, 0, 0)
    A = GeodesicWall(x1, x2, -x1)
    B = GeodesicWall(R*x1, R*x2, -R*x1)
    C = GeodesicWall(x1, x1 + x2, np.inf)
    D = GeodesicWall(x2, x1 + x2, np.inf)
    E = GeodesicWall(-x1, -x1 + x2, np.inf)
    F = GeodesicWall(-x2, x1 - x2, np.inf)

    # test plotting
    ax = axis_3d((-2, 2), (-2, 2), (0, 2))
    A.draw(ax)
    B.draw(ax)
    C.draw(ax)
    D.draw(ax)
    E.draw(ax)
    F.draw(ax)
    plt.show()

# print matrix line by line
def np_print(A) :
    for i in range(len(A)) :
        print(A[i])

# test out Mobius transformations
def test_Mobius() :
    # start with quaternion at j
    y = Quaternion(0, 0, 1, 0)
    print('y is: ' + str(y))
    print('')

    # matrices in N shift the point parallel to x1x2-plane
    n = np.array([[1, 1 + 1j], [0, 1]])
    p = mobiusTransform(n, y)
    print('n:')
    np_print(n)
    print('n(y) = ' + str(p))
    print('')

    # matrices in A move the point up the y-axis at unit speed
    a = np.array([[np.exp(2), 0], [0, np.exp(-2)]])
    p = mobiusTransform(a, y)
    print('a:')
    np_print(a)
    print('a(y) = ' + str(p))
    print('')

    # matrices in M and SO(2) fix the point
    m = np.array([[np.exp(1j*2.3), 0], [0, np.exp(-1j*2.3)]])
    p = mobiusTransform(m, y)
    print('m:')
    np_print(m)
    print('m(y) = ' + str(p))
    print('')

    o = np.array([[np.cos(-1.7), -np.sin(-1.7)], [np.sin(-1.7), np.cos(-1.7)]])
    p = mobiusTransform(o, y)
    print('o:')
    np_print(o)
    print('o(y) = ' + str(p))
    print('')

    # products of matrices act right-to-left
    X = np.matmul(np.matmul(n, a), m)
    p = mobiusTransform(X, y)
    print('nam:')
    np_print(X)
    print('nam(y) = ' + str(p))

# test Mobius transformation on a geodesic wall
def test_Mobius_wall() :
    x1 = Quaternion(1, 0, 0, 0)
    x2 = Quaternion(0, 1, 0, 0)

    R = Quaternion(np.sqrt(2), 0, 0, 0)
    A = GeodesicWall(x1, x2, -x1)
    B = GeodesicWall(R*x1, R*x2, -R*x1)
    C = GeodesicWall(x1, x1 + x2, np.inf)
    D = GeodesicWall(x2, x1 + x2, np.inf)
    E = GeodesicWall(-x1, -x1 + x2, np.inf)
    F = GeodesicWall(-x2, x1 - x2, np.inf)

    # plot original
    ax = axis_3d((-4, 4), (-4, 4), (0, 2))
    A.draw(ax)
    B.draw(ax)
    C.draw(ax)
    D.draw(ax)
    E.draw(ax)
    F.draw(ax)

    # shift everything by 2 in x2 direction
    n = np.array([[1, 1j], [0, 1]])
    A = mobius(n, A)
    B = mobius(n, B)
    C = mobius(n, C)
    D = mobius(n, D)
    E = mobius(n, E)
    F = mobius(n, F)

    # plot shifted domain
    A.draw(ax)
    B.draw(ax)
    C.draw(ax)
    D.draw(ax)
    E.draw(ax)
    F.draw(ax)
    plt.show()

# test reflections on Quaternions
# TODO: finish this!
def test_reflections() :
    x1 = Quaternion(1, 0, 0, 0)
    x2 = Quaternion(0, 1, 0, 0)

    R = Quaternion(np.sqrt(2), 0, 0, 0)
    A = GeodesicWall(x1, x2, -x1)
    B = GeodesicWall(R*x1, R*x2, -R*x1)
    C = GeodesicWall(x1, x1 + x2, np.inf)
    D = GeodesicWall(x2, x1 + x2, np.inf)
    E = GeodesicWall(-x1, -x1 + x2, np.inf)
    F = GeodesicWall(-x2, x1 - x2, np.inf)

####################
# Apollonian Group #
####################

# get the cocluster for the Apollonian group as extended geodesic walls in H3
def get_cocluster_A() :
    z = Quaternion()
    x1 = Quaternion(1)
    x2 = Quaternion(0, 1)
    h = Quaternion(1/2)
    m1 = Quaternion(-1)

    W1 = GeodesicWall(z, x2, np.inf)
    W2 = GeodesicWall(h, h + x2, np.inf)
    W3 = GeodesicWall(h*x2, x1 + h*x2, np.inf)
    W4 = GeodesicWall(x1, x2, m1*x1)

    return (W1, W2, W3, W4)

# plot the cocluster
def plot_cocluster_A(label=False) :
    # get cocluster
    Ws = get_cocluster_A()

    # plot it!
    ax = axis_3d((-1.5, 1.5), (-1.5, 1.5), (0, 2))
    for W in Ws :
        W.draw(ax, 'b')

    # if asked to label, label walls
    if label :
        ax.text(0, -1.45, 0.1, '1')
        ax.text(0.5, -1.45, 0.1, '2')
        ax.text(1.25, 0.5, 1.75, '3')
        ax.text(-0.75, -0.65, 0, '4')

    plt.show()

# get cocluster doubled across wall 1
def get_doubled_A() :
    z = Quaternion()
    x1 = Quaternion(1)
    x2 = Quaternion(0, 1)
    h = Quaternion(1/2)
    m1 = Quaternion(-1)

    Wnew = GeodesicWall(m1*h,m1*h + x2, np.inf)
    W2 = GeodesicWall(h, h + x2, np.inf)
    W3 = GeodesicWall(h*x2, x1 + h*x2, np.inf)
    W4 = GeodesicWall(x1, x2, m1*x1)

    return (Wnew, W2, W3, W4)

# plot the doubled group
def plot_doubled_A(label=False) :
    # get cocluster
    Ws = get_doubled_A()

    # plot it!
    ax = axis_3d((-1.5, 1.5), (-1.5, 1.5), (0, 2))
    for W in Ws :
        W.draw(ax, 'b')

    # if asked to label, label walls
    if label :
        ax.text(-0.5, -1.45, 0.3, '$R_1(2)$') # TODO: rotate?
        ax.text(0.5, -1.45, 0.1, '2')
        ax.text(1.25, 0.5, 1.75, '3')
        ax.text(-0.75, -0.65, 0, '4')

    plt.show()

# get flare domain for doubled group
def get_flare_A() :
    # get walls in doubled group
    Ws = get_doubled_A()

    # need to shift the fundamental domain by T for it to end up in flare
    T = np.array([ [1, 1], [0, 1] ])

    # P is the Mobius transformation mapping to flare domain
    P = np.array([ [2*np.sqrt(5), 5 - np.sqrt(5)], [-2*np.sqrt(5), 5 + np.sqrt(5)] ])

    # apply PT to each wall
    Ws_f = []
    for W in Ws :
        Ws_f.append(mobius(P, mobius(T, W)))

    return Ws_f

# plot flare domain for doubled group
def plot_flare_A() :
    # grab domain
    Ws = get_flare_A()

    # plot it!
    ax = axis_3d((-5, 5), (-5, 5), (0, 9))
    for W in Ws :
        W.draw(ax, 'b')

    plt.show()

##################
# Pullback Tests #
##################

# pullback algorithm for the double Apollonian group (also keeps track of word
# in matrices used to map to pullback)
def pullback_A(z) :
    # here are the matrices involved in the pullback algorithm
    T = np.array([ [1, 1],
                   [0, 1]])
    Tinv = np.array([ [1, -1],
                      [0, 1]])
    A = np.array([ [1j, 1],
                   [0, -1j]])
    S = np.array([ [0, -1],
                   [1, 0]])

    # keep track of word in matrices
    word = ''

    # first, need to apply A if x2 is less than 1/2
    if z.x2 < 1/2 :
        z = mobiusTransform(A, z)
        word = 'A'

    # now repeat strings of T's then S until done
    while 1 :
        # use T's to get x1 between -1/2 and 1/2
        numTs = 0
        while z.x1 < -1/2 :
            z = mobiusTransform(T, z)
            numTs += 1
        while z.x1 >= 1/2 :
            z = mobiusTransform(Tinv, z)
            numTs -= 1

        # add the T-string to word
        if numTs < 0 :
            word = 'T^(' + str(numTs) + ')' + word
        elif numTs > 0 :
            word = 'T^' + str(numTs) + word

        # if |z| >= 1, done; otherwise, apply S
        if z.norm() >= 1 :
            return z, word
        else :
            z = mobiusTransform(S, z)
            word = 'S' + word

# map from original space to flare domain
def to_flare(z) :
    # P is the Mobius transformation mapping to flare domain
    P = np.array([ [2*np.sqrt(5), 5 - np.sqrt(5)], [-2*np.sqrt(5), 5 + np.sqrt(5)] ])

    # need to shift point by T for it to end up in flare
    return mobiusTransform(P, Quaternion(z.x1 + 1, z.x2, z.y))

# map from flare domain to original space
def to_orig(z) :
    # B is the Mobius transformation mapping back from flare domain
    B = np.array([ [-5 - np.sqrt(5), 5 - np.sqrt(5)], [-2*np.sqrt(5), -2*np.sqrt(5)] ])

    # shift back by T^(-1) at end
    z_temp = mobiusTransform(B, z)
    return Quaternion(z_temp.x1 - 1, z_temp.x2, z_temp.y)

# playing with test points and pullbacks
def test_test_points(N, thet, phi) :
    # scaling parameter for Apollonian packing
    kappa = (3 + np.sqrt(5))/(3 - np.sqrt(5));

    # form test points in flare domain
    rhos = np.exp(np.arange(1, N + 1)/(2*N)*np.log(kappa))
    z_flares = []
    for rho in rhos :
        z_flares.append(Quaternion(rho*np.sin(phi)*np.cos(thet),
                                   rho*np.sin(phi)*np.sin(thet),
                                   rho*np.cos(phi)))

    # get points in original domain plus pullbacks (and look at word for each
    # pullback)
    zs = []
    z_pbs = []
    for z in z_flares :
        z_orig = to_orig(z)
        z_pb, word = pullback_A(z_orig)
        print(word)

        zs.append(z_orig)
        z_pbs.append(z_pb)

    # plot some points
    # Ws = get_doubled_A()
    # ax = axis_3d((-3, 3), (-3, 3), (0, 5))
    # for W in Ws :
    #     W.draw(ax, 'b')
    # plot_points(ax, zs, 'o', color='k')
    # plot_points(ax, z_pbs, 'o', color='g')
    #
    # plt.show()

# create list of thetas and phis to play with
def many_test_points() :
    N = 15
    thets = np.linspace(0, 2*np.pi - 0.05, 20)
    phis = np.linspace(np.pi/4, np.pi/2 - 0.01, 3)

    for thet in thets :
        for phi in phis :
            print('thet =', thet, 'phi =', phi)
            test_test_points(N, thet, phi)
            print('')

if __name__ == '__main__' :
    plot_flare_A()
