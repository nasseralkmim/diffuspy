import matplotlib.pyplot as plt
from diffuspy.postprocess import draw
import matplotlib.animation as animation
import numpy as np


def show():
    plt.show()


def model(model, name=None, color='k', dpi=100, ele=False, ele_label=False,
          surf_label=False, nodes_label=False, edges_label=False, surf=False):
    """Plot the  model geometry

    """
    fig = plt.figure(name, dpi=dpi)
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$y$')
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
        # # Change the colorbar range
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=plt.Normalize(vmin=np.amin(U),
                                                      vmax=np.amax(U)))
        # # fake up the array of the scalar mappable. Urgh...
        sm._A = []
        cbar = plt.colorbar(sm)
        cbar.set_label(r'Temperature $^{\circ}C$')

    return cs


def contour_animation(model, T, t_int, dt, name='Temperature.gif',
                      bitrate=300, lev=10, interval=100,
                      time_text_color='black'):
    """Plot the animation for the temperature evolution countor plot

    """
    print('Initializing plotter...')
    N = int(t_int/dt)

    frm = []

    fig, ax = plt.subplots()

    vmin, vmax = np.amin(T), np.amax(T)

    for n in range(N+1):
        t = n * dt
        im = draw.tricontourf(model, T[:, n], ax, cmap='hot', lev=lev,
                              vmin=vmin, vmax=vmax, cl=False)

        te = ax.text(0, 1, "Time (h): "+str(round(t/(60*60), 2)),
                     ha='left', va='top',
                     transform=ax.transAxes, color=time_text_color)
        frm.append(im.collections + [te])

    # Change the colorbar range
    sm = plt.cm.ScalarMappable(cmap='hot',
                               norm=plt.Normalize(vmin=vmin, vmax=vmax))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    cbar = plt.colorbar(sm)
    cbar.set_label(r'Temperature $^{\circ} C$')

    ani = animation.ArtistAnimation(fig, frm, interval=interval, blit=True)
    ani.save(name, writer='imagemagick', bitrate=bitrate)
    print('Plotting completed!')


def function(f, interval, line, xlabel, ylabel, **kwargs):
    """Plots the function in the interval

    """
    fig, ax = plt.subplots()

    x = np.linspace(0, interval)

    ax.plot(x/3600, f(1, 1, x)[line], **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.tight_layout()
