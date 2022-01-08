from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge

class Geodesic :
    """A single geodesic in the upper half plane model of hyperbolic space."""

    # initialize a Geodesic by specifying the two end points on the extended reals
    def __init__(self, endpoint1, endpoint2) :
        # verify input
        if endpoint1 == endpoint2 :
            raise ValueError('Endpoints of a geodesic must be distinct.')
        if endpoint1 != np.inf :
            if endpoint1.imag != 0 :
                raise ValueError('Endpoints must be extended reals.')
        if endpoint2 != np.inf :
            if endpoint2.imag != 0 :
                raise ValueError('Endpoints must be extended reals.')

        # initialize variables
        self.ep1 = endpoint1
        self.ep2 = endpoint2

    # draw this geodesic on an Axes object from matplotlib.artist.Artist
    def draw(self, ax, hue='k', lstyle='-', angle=180) :
        yMax = plt.ylim()[1]

        if self.ep1 == np.inf :
            ax.plot( [self.ep2, self.ep2], [0, yMax], color=hue, ls=lstyle )
        elif self.ep2 == np.inf :
            ax.plot( [self.ep1, self.ep1], [0, yMax], color=hue, ls=lstyle )
        else :
            diam = abs(self.ep1 - self.ep2)
            wedge = Wedge( ((self.ep1 + self.ep2)/2.0, 0), diam/2.0, 0, angle, color=hue, ls=lstyle, fill=0)
            ax.add_patch(wedge)

    # return string with the endpoints of the geodesic
    def info(self) :
        return 'Geodesic: ' + str(self.ep1) + ', ' + str(self.ep2)

class FundamentalDomain :
    """A fundamental domain for a Fuchsian group in the upper half plane model of hyperbolic space.
       This is described by some positive number of geodesics."""

    # initialize a FundamentalDomain by specifying a list of geodesics
    def __init__(self, gs, refPnt=None) :
        self.geodesics = gs
        self.referencePoint = refPnt

    # draw this fundamental domain on an Axes object
    def draw(self, ax, hue='k', lstyle='-', fill=0, angle=180) :
        for g in self.geodesics :
            g.draw(ax, hue, lstyle, angle)

        # fill in fundamental domain if fill is 1
        if fill :
            # can't fill without a reference point
            if self.referencePoint == None :
                raise ValueError("Can't fill in fundamental domain without a reference point.")

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
                # we're either to the left or right of a vertical geodesic
                if g.ep1 == np.inf or g.ep2 == np.inf :
                    xVal = min(g.ep1, g.ep2)
                    if xRef >= xVal :
                        Zs.append(np.greater_equal(X, xVal*N))
                    else :
                        Zs.append(np.less(X, xVal*N))

                # we're either inside or outside of a semi-circular geodesic
                else :
                    xCent = (g.ep1 + g.ep2)/2
                    rsq = (g.ep2 - xCent)**2
                    if (xRef - xCent)**2 + yRef**2 >= rsq :
                        Zs.append( np.greater_equal((X - xCent)**2 + Y**2, rsq*N) )
                    else :
                        Zs.append( np.less((X - xCent)**2 + Y**2, rsq*N) )

            Z = np.greater_equal(Y, 0*N)
            for i in range(0, len(Zs)) :
                Z = Z*Zs[i]

            # use imshow to draw colors
            ax.imshow(Z, cmap='gray_r', alpha=0.25, extent=[xMin,xMax,yMax,yMin])

    # print all geodesic info for this fundamental domain
    def print_info(self) :
        print('Fundamental Domain')
        for g in self.geodesics :
            print('    ' + g.info())

##########################
# Mobius Transformations #
##########################

# given invertible 2x2 matrix A, returns A(z) where action is by Mobius transformation
def mobiusTransform(A, z) :
    if z == np.inf :
        if A[1, 0] == 0.0 :
            return np.inf
        return A[0, 0] / A[1, 0]
    elif A[1, 0] * z + A[1, 1] == 0.0 :
        return np.inf
    return (A[0, 0] * z + A[0, 1]) / (A[1, 0] * z + A[1, 1])

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

# sets up figure of specified dimensions with real line drawn
def setupFig(xMin, xMax, yMin, yMax, hue='k') :
    fig = plt.figure()
    ax = plt.axes( xlim=(xMin, xMax), ylim=(yMin, yMax), aspect='equal')
    ax.plot([xMin, xMax], [0, 0], color=hue)
    return ax

############
# Examples #
############

# creates a FundamentalDomain object for the reflection group whose double
# cover (when doubled across the imaginary axis) is the classic fundamental
# domain for SL(2, Z)
def refl_group_FD() :
    gs = [Geodesic(0, np.inf), Geodesic(1/2, np.inf), Geodesic(-1, 1)]
    return FundamentalDomain(gs, 1/4 + 2j)

# creates a FundamentalDomain object for a Hecke triangle group
#   - we assume the group is generated by z -> z+1 and z -> -r^2/z
#   - this assumption is conjugate to the case of z -> z+mu, z -> -1/z^2
def heckeFD(r) :
    gs = [Geodesic(0, np.inf), Geodesic(1, np.inf), Geodesic(-r, r), Geodesic(1-r, 1+r)]
    return FundamentalDomain(gs, 1/2 + 1j)

# creates a FundamentalDomain object for a Hecke triangle group
#   - we assume the group is generated by z -> z+2mu and z -> -1/z
#   - this assumption is conjugate to the case of z -> z+1, z -> -R^2/z^2
def heckeGroupFD(mu) :
    gs = [Geodesic(-mu, np.inf), Geodesic(mu, np.inf), Geodesic(-1, 1)]
    return FundamentalDomain(gs, 2j)

# plots flare version of the Hecke group as obtained from heckeFD
def plot_Hecke_flare(r) :
    # get original fundamental domain
    fd = heckeFD(r)

    # z1 and z2 are endpoints of the axis of our hyperbolic element
    z1 = (1 - np.sqrt(1 - 4*r**2))/2
    z2 = (1 + np.sqrt(1 - 4*r**2))/2

    # set up matrix to map to flare
    c = (r - z2)/(r - z1)
    U = np.array([[c, -c*z1], [1, -z2]])

    # now actually map to flare!
    fd = mobius(U, fd)

    # plot
    kappa = (z2/r)**2

    ax = setupFig(-7, 7, 0, 9)
    fd.draw(ax, fill=1)

    # labels
    plt.xticks([1, kappa], ['$1$', '$\kappa$'])
    plt.yticks([])

    # plot axis at angle alpha as dotted line
    R = 15
    alpha = np.pi/6
    ax.plot([0, R*np.cos(alpha)], [0, R*np.sin(alpha)], 'k--')

    plt.show()

# creates and draws a FundamentalDomain object for the Apollonian cocluster
# intersected with the x1x2-plane and doubled across the x2-axis
def draw_apollonian_FD() :
    # most lines
    gs = [Geodesic(0.5, np.inf),
          Geodesic(-0.5, np.inf),
          Geodesic(-1, 1)]
    fd = FundamentalDomain(gs, 2j)

    # set up figure and draw geodesics
    ax = setupFig(-2, 2, -2, 2)
    fd.draw(ax, hue='b', fill=1, angle=360)

    # draw missing horizontal line and rest of vertical lines
    ax.plot( [-2, 2], [0.5, 0.5], color='b')
    ax.plot( [-0.5, -0.5], [-2, 2], color='b')
    ax.plot( [0.5, 0.5], [-2, 2], color='b')
    ax.plot( [0, 0], [-2, 2], 'b--', alpha=0.5)

    # label everything
    plt.text(1.25, -0.1, 'v')
    plt.text(1.25, 0.55, '3')
    plt.text(-0.75, 0.75, '4')
    plt.text(0.55, -1.25, '2')
    plt.text(-0.1, -1.25, '1')

    plt.show()

# draws a flare example, before being mapped to flare domain
def pre_flare() :
    g1 = Geodesic(-1, -1/4)
    g2 = Geodesic(1/4, 1)
    g = Geodesic(-1/2, 1/2)

    ax = setupFig(-1.5, 1.5, 0, 1)

    g1.draw(ax)
    g2.draw(ax)
    g.draw(ax, 'k', ':')

    plt.xticks([-1/2, -1/4, 1/2, 1], ['$z_1$', '$t$', '$z_2$', '$s$'])
    plt.yticks([])

    plt.show()

# draws a flare example, after being mapped to flare domain
def post_flare() :
    g1 = Geodesic(-1, 1)
    g2 = Geodesic(-2, 2)
    g = Geodesic(0, np.inf)

    f = FundamentalDomain([g1, g2], 1.5j)

    ax = setupFig(-3, 3, 0, 3)

    f.draw(ax, fill=1)
    #g.draw(ax, 'k', ':')

    #plt.xticks([0, 1], ['$U(z_1) = 0$', '$U(t) = 1$'])
    plt.xticks([1, 2], ['$1$', '$\kappa$'])
    plt.yticks([])

    plt.show()

if __name__ == '__main__' :
    post_flare()
