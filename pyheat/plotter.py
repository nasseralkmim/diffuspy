import matplotlib.pyplot as plt
from elastopy import draw


def show():
    plt.show()


def model(mesh, name=None, color='k', dpi=100, ele=False, ele_label=False,
          surf_label=False, nodes_label=False, edges_label=False):
    """Plot the  model geometry

    """
    fig = plt.figure(name, dpi=dpi)
    ax = fig.add_subplot(1, 1, 1)
    ax.set_xlabel(r'x')
    ax.set_ylabel(r'y')
    ax.set_aspect('equal')

    draw.domain(mesh, ax, color=color)

    if ele is True:
        draw.elements(mesh, ax, color=color)

    if ele_label is True:
        draw.elements_label(mesh, ax)

    if surf_label is True:
        draw.surface_label(mesh, ax)

    if nodes_label is True:
        draw.nodes_label(mesh, ax)

    if edges_label is True:
        draw.edges_label(mesh, ax)

    return None


def contour(mesh, U, cmap='hot', lev=10):
    """Plot stress with nodal stresses

    """
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_aspect('equal')
    ax.set_xlabel(r'x')
    ax.set_ylabel(r'y')

    ax.set_title(r'Temperature (C)')
    draw.tricontourf(mesh, U, ax, cmap=cmap, lev=lev)

    return None
