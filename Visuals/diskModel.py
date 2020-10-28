from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge
from matplotlib.patches import Arc

###########################################
# Geometric helper functions
###########################################

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

# computes Cayley transform of z
def cayley_transform(z, r) :
    if z == np.inf :
        return 1
    return (z - r*1j) / (z + r*1j)

###########################################
# Geodesic and FundamentalDomain classes
###########################################

class Geodesic :
    """A single geodesic in the disk model of hyperbolic space."""

    # initialize a Geodesic by specifying the two end points on the unit circle
    def __init__(self, endpoint1, endpoint2) :
        # verify input
        if endpoint1 == endpoint2 :
            raise ValueError('Endpoints of a geodesic must be distinct.')
        if abs(abs(endpoint1) - 1) > 0.001 :
            raise ValueError('Endpoints must on the unit circle.')
        if abs(abs(endpoint2) - 1) > 0.001 :
            raise ValueError('Endpoints must on the unit circle.')

        # initialize variables
        self.ep1 = endpoint1
        self.ep2 = endpoint2

    # draw this geodesic on an Axes object from matplotlib.artist.Artist
    def draw(self, ax, hue='k', lstyle='-') :
        # get angles of endpoints
        theta1 = rectToPol(self.ep1.real, self.ep1.imag)[1]
        theta2 = rectToPol(self.ep2.real, self.ep2.imag)[1]
        if theta1 == theta2 :
            raise Exception('theta1 cannot equal theta2')

        # find acute angle between thetas
        ang = np.abs(theta1 - theta2)
        wasObtuse = 0
        if ang == np.pi :
            # plot line between them
            x1, y1 = polToRect(1, theta1)
            x2, y2 = polToRect(1, theta2)
            ax.plot( [x1, x2], [y1, y2], color=hue, ls=lstyle)
            return
        elif ang > np.pi :
            ang = 2*np.pi - ang
            wasObtuse = 1

        # Assume we rotate to be symmetric across x-axis and in 1st/4th quadrant
        theta = ang / 2

        # get coordinates of center
        rC = 1 / np.cos(theta)
        thetaC = (theta1 + theta2) / 2
        if wasObtuse :
            thetaC += np.pi
        xC, yC = polToRect(rC, thetaC)

        # get diameter
        xP, yP = polToRect(1, theta)
        diam = 2*np.sqrt( (xP - rC)**2 + yP**2 )

        # plot arc
        ang1 = np.pi/2 + theta + thetaC
        ang2 = 3*np.pi/2 - theta + thetaC
        arc = Arc( (xC, yC), diam, diam, theta1=180*ang1/np.pi, theta2=180*ang2/np.pi, color=hue, ls=lstyle)
        ax.add_patch(arc)

class FundamentalDomain :
    """A fundamental domain for a Fuchsian group in the upper half plane model of hyperbolic space.
       This is described by some positive number of geodesics."""

    # initialize a FundamentalDomain by specifying a list of geodesics
    def __init__(self, gs, refPnt=None) :
        self.geodesics = gs
        self.referencePoint = refPnt

    # draw this fundamental domain on an Axes object
    def draw(self, ax, hue='k', lstyle='-') :
        for g in self.geodesics :
            g.draw(ax, hue, lstyle)

        # fill in fundamental domain if reference point is given
        if self.referencePoint != None :
            # get axis bounds and real/imag parts of reference point
            xMin, xMax = ax.get_xlim()
            yMin, yMax = ax.get_ylim()
            xRef = (self.referencePoint).real
            yRef = (self.referencePoint).imag

            # set up boolean matrix for whether points are in fundamental domain
            delta = min( (xMax-xMin)/1000, (yMax-yMin)/1000 )
            x = np.arange(xMin, xMax, delta)
            y = np.arange(yMin, yMax, delta)
            X, Y = np.meshgrid(x, y)
            N = np.ones(X.shape)
            Zs = []
            for g in self.geodesics :
                # get angles of endpoints
                theta1 = rectToPol(g.ep1.real, g.ep1.imag)[1]
                theta2 = rectToPol(g.ep2.real, g.ep2.imag)[1]

                # need angle to tell if geodesic is line or circle
                ang = np.abs(theta1 - theta2)

                # we're on either side of a line
                if ang == np.pi :
                    # rotate by -theta1 so it's a horizontal line
                    if -np.sin(theta1)*xRef + np.cos(theta1)*yRef < 0 :
                        Zs.append(np.less(-np.sin(theta1)*X + np.cos(theta1)*Y, 0))
                    else :
                        Zs.append(np.greater_equal(-np.sin(theta1)*X + np.cos(theta1)*Y, 0))

                # we're either inside or outside of a semi-circular geodesic
                else :
                    wasObtuse = 0
                    if ang > np.pi :
                        ang = 2*np.pi - ang
                        wasObtuse = 1

                    # Assume we rotate to be symmetric across x-axis and in 1st/4th quadrant
                    theta = ang / 2

                    # get coordinates of center
                    rC = 1 / np.cos(theta)
                    thetaC = (theta1 + theta2) / 2
                    if wasObtuse :
                        thetaC += np.pi
                    xC, yC = polToRect(rC, thetaC)

                    # get radius squared
                    xP, yP = polToRect(1, theta)
                    rsq = (xP - rC)**2 + yP**2

                    # fill points
                    if (xRef - xC)**2 + (yRef - yC)**2 >= rsq :
                        Zs.append( np.greater_equal((X - xC)**2 + (Y - yC)**2, rsq*N) )
                    else :
                        Zs.append( np.less((X - xC)**2 + (Y - yC)**2, rsq*N) )

            Z = np.less_equal(X**2 + Y**2, N)
            for i in range(0, len(Zs)) :
                Z = Z*Zs[i]

            # use imshow to draw colors
            ax.imshow(Z, cmap='gray_r', alpha=0.25, extent=[xMin,xMax,yMax,yMin])

###########################################
# Testing and related functions
###########################################

# sets up figure with unit circle drawn
def setupFig(hue='k') :
    fig = plt.figure()
    ax = plt.axes( xlim=(-1.1, 1.1), ylim=(-1.1, 1.1), aspect='equal')
    unitCirc = Arc( (0, 0), 2, 2, theta1=0, theta2=360, color='k')
    ax.add_patch(unitCirc)
    return ax

# creates a FundamentalDomain object for a Hecke triangle group
#   - we assume the group is generated by z -> z+1 and z -> -r^2/z
#   - this assumption is conjugate to the case of z -> z+mu, z -> -1/z^2
def heckeFD(r, type=1) :
    if type == 1 :
        g1 = Geodesic(cayley_transform(-1/2, r), cayley_transform(np.inf,r ))
        g2 = Geodesic(cayley_transform(1/2, r), cayley_transform(np.inf, r))
        g3 = Geodesic(cayley_transform(-r, r), cayley_transform(r, r))
        gs = [g1, g2, g3]
        return FundamentalDomain(gs, cayley_transform(2j, r))
    else :
        g1 = Geodesic(cayley_transform(0, r), cayley_transform(np.inf, r))
        g2 = Geodesic(cayley_transform(1, r), cayley_transform(np.inf, r))
        g3 = Geodesic(cayley_transform(-r, r), cayley_transform(r, r))
        g4 = Geodesic(cayley_transform(1-r, r), cayley_transform(1+r, r))
        gs = [g1, g2, g3, g4]
        return FundamentalDomain(gs, cayley_transform(1/2 + 1j, r))

# creates a FundamentalDomain object for a symmetric Schottky group
#   - the group is generated by reflections through N circles equally
#     spaced around the circle, cutting out an angle of theta
#   - note that theta must be less than or equal to 2pi/N for the circles
#     to not overlap
def schottky(N, theta) :
    if theta > 2*np.pi/N :
        raise ValueError("theta must be <= 2pi/N")

    # get angles for centers of bounding circles
    center_angles = [0]
    for i in range(1, N) :
        center_angles.append(i*2*np.pi/N)

    # now get list of geodesics
    gs = []
    for c in center_angles :
        gs.append(Geodesic(np.exp(1j*(c - theta/2)), np.exp(1j*(c + theta/2))))

    return FundamentalDomain(gs, 0)

def test1() :
    fd = heckeFD(7/20, 2)

    ax = setupFig()
    fd.draw(ax)

    # remove axes and show figure
    plt.axis("off")
    plt.show()

def test2() :
    ax = setupFig()
    N = 15
    for i in range(N) :
        z = cayley_transform(i/(2*N) + 0.32j, r)
        plt.plot([z.real], [z.imag], 'k.')

    # remove axes and show figure
    plt.axis("off")
    plt.show()

def make_schottky_pic() :
    # set up and draw Schottky domain
    fd = schottky(3, np.pi/2)
    ax = setupFig()
    fd.draw(ax)

    # add wedge and theta angle
    wedge = Wedge( (0, 0), 1, theta1=-45, theta2=45, color='k', fill=0)
    ax.add_patch(wedge)
    plt.text(0.05, -0.03, "$\\theta$")

    # remove axes and show figure
    plt.axis("off")
    plt.show()

if __name__ == '__main__' :
    make_schottky_pic()
