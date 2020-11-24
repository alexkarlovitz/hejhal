from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import poincareModel as pm
import diskModel as dm
from schottkyGroups import cayley_transform, inv_cayley_transform, \
                           reflect, get_UHP_circles, get_FD, get_axes

# plot Poincare model fundamental domain in the disk model on axis ax
def plot_in_disk(fd, ax) :
    # create list of corresponding geodesics in disk model
    gs = []
    for g in fd.geodesics :
        gs.append(dm.Geodesic(cayley_transform(g.ep1), cayley_transform(g.ep2)))

    # create new fundamental domain object and plot
    refPoint = cayley_transform(fd.referencePoint)
    fd_new = dm.FundamentalDomain(gs, refPoint)
    fd_new.draw(ax)

# reflects each geodesic in fd by Rs[i] for i in ts
def map_fd(fd, ts, Rs) :
    # get copy of geodesics
    gs = []
    for g in fd.geodesics :
        gs.append(pm.Geodesic(g.ep1, g.ep2))

    # map each geodesic and reference point by each reflection
    ref_point = fd.referencePoint
    for i in ts :
        for j in range(len(gs)) :
            gs[j] = pm.Geodesic(reflect(gs[j].ep1, Rs[i][0], Rs[i][1]), \
                                reflect(gs[j].ep2, Rs[i][0], Rs[i][1]))

        # also map reference point
        ref_point = reflect(ref_point, Rs[i][0], Rs[i][1])

    return pm.FundamentalDomain(gs, ref_point)

def test_mapped_fd() :
    circles = get_UHP_circles(3, np.pi/6)
    fd = get_FD(circles)

    fd_new = map_fd(fd, [0, 1, 2, 1, 2], circles)

    ax = dm.setupFig()
    plot_in_disk(fd_new, ax)

    _, _, g = get_axes(np.pi/6)
    gd = dm.Geodesic(cayley_transform(g.ep1), cayley_transform(g.ep2))
    gd.draw(ax)

    plt.show()

if __name__ == '__main__' :
    test_mapped_fd()
