from __future__ import division
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

    # draw this geodesic wall on a given matplotlib axis
    #   see http://web.archive.org/web/20161011113446/http://www.abecedarical.com/zenosamples/zs_circle3pts.html
    def draw(self, ax) :
        # check if one of the points is inf
        if self.q1 == np.inf or self.q2 == np.inf or self.q3 == np.inf :
            # find non-infinite points
            p1, p2 = self.q1, self.q2
            if p1 == np.inf :
                p1 = self.q3
            elif p2 == np.inf :
                p2 = self.q3

            # draw the half plane
            self.draw_plane(ax, p1.x1, p1.x2, p2.x1, p2.x2)

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
                self.draw_plane(ax, x1, y1, x2, y2)

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

                # draw half sphere
                self.draw_sphere(ax, x, y, r)

    # given two points on x1x2-plane, draw vertical half plane through them
    def draw_plane(self, ax, x1, y1, x2, y2) :
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
            ax.plot_surface(X, Y, Z, color='k', alpha=0.2)

        # check if plane is parallel to x2 axis
        elif x1 == x2 :
            # make mesh in y and z coordinates
            Y, Z = np.mgrid[y_min:y_max:2j, z_min:z_max:2j]

            # now plot
            X = np.zeros_like(Y) + x1
            ax.plot_surface(X, Y, Z, color='k', alpha=0.2)

        # otherwise, treat y as a function of x
        else :
            # make mesh in x and z coordinates
            X, Z = np.mgrid[x_min:x_max:2j, z_min:z_max:2j]

            # compute y and plot
            m = (y2 - y1)/(x2 - x1)
            Y = m*(X - x1) + y1
            ax.plot_surface(X, Y, Z, color='k', alpha=0.2)

    # given center (on x1x2-plane) and radius, draw half sphere
    def draw_sphere(self, ax, x, y, r) :
        # form meshgrid for unit half sphere, then scale and shift
        u, v = np.mgrid[0:2*np.pi:51j, 0:np.pi/2:51j]
        X = r*np.cos(u)*np.sin(v) + x
        Y = r*np.sin(u)*np.sin(v) + y
        Z = r*np.cos(v)

        # plot it!
        ax.plot_surface(X, Y, Z, color='k', alpha=0.2)

    # return string with the endpoints of the geodesic
    def __str__(self) :
        return str(self.q1) + ', ' + str(self.q2) + ', ' + str(self.q3)

#####################
# Mobius Transforms #
#####################

# TODO: double check this func, then update one below

# given invertible 2x2 matrix A, returns A(z) where action is by Mobius transformation
#   - A has entries in C
#   - z is a Quaternion or np.inf
def mobiusTransform(A, z) :
    # convert entries in A to Quaternions
    B = [[Quaternion(A[0, 0].real, A[0, 0].imag), Quaternion(A[0, 1].real, A[0, 1].imag)],
         [Quaternion(A[1, 0].real, A[1, 0].imag), Quaternion(A[1, 1].real, A[1, 1].imag)]]

    if z == np.inf :
        if B[1][0] == 0.0 :
            return np.inf
        return B[0][0] / B[1][0]
    elif B[1][0] * z + B[1][1] == 0.0 :
        return np.inf
    return (B[0][0] * z + B[0][1]) / (B[1][0] * z + B[1][1])

# given invertible 2x2 matrix A, returns mobius transformation A of a geodesic or fund domain
def mobius(A, x) :
    # if x is a geodesic, get new geodesic where we applied A to endpoints
    if isinstance(x, Geodesic) :
        e1 = mobiusTransform(A, x.ep1)
        e2 = mobiusTransform(A, x.ep2)
        return Geodesic(e1, e2)

    # if x is a fundamental domain, get new fundamental domain by applying A to all geodesics
    elif isinstance(x, FundamentalDomain) :
        gsNew = []
        for g in x.geodesics :
            gsNew.append( mobius(A, g) )
        refNew = None
        if x.referencePoint != None :
            refNew = mobiusTransform(A, x.referencePoint)
        return FundamentalDomain(gsNew, refNew)

    else :
        raise ValueError('Function "mobius" expects second input to be a Geodesic or FundamentalDomain object.')

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
        print A[i]

# test out Mobius transformations
def test_Mobius() :
    # start with quaternion at j
    y = Quaternion(0, 0, 1, 0)
    print 'y is: ' + str(y)
    print ''

    # matrices in N shift the point parallel to x1x2-plane
    n = np.array([[1, 1 + 1j], [0, 1]])
    p = mobiusTransform(n, y)
    print 'n:'
    np_print(n)
    print 'n(y) = ' + str(p)
    print ''

    # matrices in A move the point up the y-axis at unit speed
    a = np.array([[np.exp(2), 0], [0, np.exp(-2)]])
    p = mobiusTransform(a, y)
    print 'a:'
    np_print(a)
    print 'a(y) = ' + str(p)
    print ''

    # matrices in M and SO(2) fix the point
    m = np.array([[np.exp(1j*2.3), 0], [0, np.exp(-1j*2.3)]])
    p = mobiusTransform(m, y)
    print 'm:'
    np_print(m)
    print 'm(y) = ' + str(p)
    print ''

    o = np.array([[np.cos(-1.7), -np.sin(-1.7)], [np.sin(-1.7), np.cos(-1.7)]])
    p = mobiusTransform(o, y)
    print 'o:'
    np_print(o)
    print 'o(y) = ' + str(p)
    print ''

    # products of matrices act right-to-left
    X = np.matmul(np.matmul(n, a), m)
    p = mobiusTransform(X, y)
    print 'nam:'
    np_print(X)
    print 'nam(y) = ' + str(p)

if __name__ == '__main__' :
    test_Mobius()
