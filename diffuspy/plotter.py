import matplotlib.pyplot as plt
from diffuspy import draw
import matplotlib.animation as animation


def show():
    plt.show()


def model(model, name=None, color='k', dpi=100, ele=False, ele_label=False,
          surf_label=False, nodes_label=False, edges_label=False, surf=False):
    """Plot the  model geometry

    """
    fig = plt.figure(name, dpi=dpi)
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.set_xlabel(r'x')
    ax.set_ylabel(r'y')
    ax.set_aspect('equal')

    draw.domain(model, ax, color=color)

    if ele is True:
        draw.elements(model, ax, color=color)

    if ele_label is True:
        draw.elements_label(model, ax)

    if surf_label is True:
        draw.surface_label(model, ax)

    if nodes_label is True:
        draw.nodes_label(model, ax)

    if edges_label is True:
        draw.edges_label(model, ax)

    if surf is True:
        draw.surface(model, ax)

    return None


def contour(model, U, cmap='hot', lev=10, name=None, contour_label=True,
            cbar=True, figure=True, vmin=None, vmax=None, title=''):
    """Plot stress with nodal stresses

    """
    if figure is True:
        fig = plt.figure(name)
        ax = fig.add_axes([.1, .1, .8, .8])
        ax.set_aspect('equal')
        ax.axis('off')
    else:
        ax = plt.gca()

    ax.set_title(title)

    cs = draw.tricontourf(model, U, ax, cmap=cmap, lev=lev,
                          cl=contour_label, vmin=vmin, vmax=vmax)

    if cbar is True:
        # cbar = fig.colorbar(cs)
        # # Change the colorbar range
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=plt.Normalize(vmin=vmin, vmax=vmax))
        # # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        cbar = plt.colorbar(sm)
        cbar.set_label(r'Temperature $^{\circ}C$')
    return cs


def anime(frames, dt):
    """Plot animation with images frames

    """
    fig = plt.gcf()
    ani = animation.ArtistAnimation(fig, frames, interval=dt, blit=True)
    return ani
