import matplotlib.pyplot as plt
from diffuspy.postprocess import draw
import matplotlib.animation as animation
import numpy as np


def show():
    plt.show()


def model(model, name=None, color='k', dpi=100, ele=False, ele_label=False,
          surf_label=False, nodes_label=False, edges_label=False, surf=False,
          font_size=12,
          ax=None):
    """Plot the  model geometry

    """
    if ax is None:
        fig = plt.figure(name, dpi=dpi)
        ax = fig.add_axes([.1, .1, .8, .8])

    ax.set_xlabel(r'$x (m)$')
    ax.set_ylabel(r'$y (m)$')
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


def contour(model, U, cmap='hot', lev=10, name=None, ax=None,
            contour_label=True,
            cbar=True, vmin=None, vmax=None, title='',
            font_size=12, axis='off'):
    """Plot stress with nodal stresses

    """
    print('Initializing Plotter...', end='')
    if ax is None:
        fig = plt.figure(name)
        ax = fig.add_axes([.1, .1, .8, .8])
        ax.set_aspect('equal')
        ax.axis(axis)
    else:
        ax = ax
        ax.set_aspect('equal')
        ax.axis(axis)

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

    print('Plotting done!')

    return cs


def contour_animation(model, T, t_int, dt, name='Temperature.gif',
                      bitrate=300, lev=10, interval=100,
                      time_text_color='black', font_size=12,
                      time_scale='hour'):
    """Plot the animation for the temperature evolution countor plot

    """
    print('Plotting animation...', end='')
    N = int(t_int/dt)

    if time_scale == 'day':
        time_factor = 24*60*60
    else:
        time_factor = 60*60

    frm = []

    fig, ax = plt.subplots()
    ax.set_xlabel(r'$x$', size=font_size)
    ax.set_ylabel(r'$y$', size=font_size)
    ax.set_aspect('equal')

    vmin, vmax = np.amin(T), np.amax(T)

    for n in range(N+1):
        t = n * dt
        im = draw.tricontourf(model, T[:, n], ax, cmap='hot', lev=lev,
                              vmin=vmin, vmax=vmax, cl=False)

        te = ax.text(0, 1,
                     "Time (" + time_scale + "): " +
                    str(round(t/(time_factor), 2)),
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
    print('Completed!')


def function(f, interval, line, xlabel, ylabel, **kwargs):
    """Plots the function in the interval

    """
    fig, ax = plt.subplots()

    x = np.linspace(0, interval)

    ax.plot(x/3600, f(1, 1, x)[line], **kwargs)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.tight_layout()


def array(T, node, interval, dt, font_size=12, **kwargs):
    """Plot the array at a specific point in time

    """
    fig, ax = plt.subplots()

    t = np.linspace(0, interval/3600, interval/dt+1)

    ax.plot(t, T[node, :], **kwargs)
    ax.set_xlabel('Time, h', size=font_size)
    ax.set_ylabel(r'Temperature $^{\circ}C$', size=font_size)
    plt.tight_layout()


def solution_at_time(model, T, time, t_int, dt, time_scale='day',
                     ax=None, y=0, **plot_kwargs):
    """Plot solution at a specific time

    """
    print('Plotting solution at time {} through the x-axis '.format(time), end='')
    if ax is None:
        fig, ax = plt.subplots()

    if time_scale == 'day':
        # (t_int/dt) is the number of steps
        # (t_int/60*60*24) is the time interval in days
        time_index = int((t_int/dt)/(t_int/(60*60*24))*time)
    elif time_scale == 'hour':
        time_index = int((dt/t_int)/(t_int/60*60)*time)
    else:
        print('Time scale should be hour or day!')

    nodes = np.where(model.XYZ[:, 1] == y)[0]
    x = model.XYZ[nodes, 0]
    data = np.array(list(zip(x, T[nodes, time_index])))

    sorted_data = data[np.argsort(data[:, 0])]

    ax.plot(sorted_data[:, 0], sorted_data[:, 1], **plot_kwargs)
    ax.set_xlabel(r'$(m)$')
    ax.set_ylabel(r'Temperature $^{\circ}C$')
    plt.legend()
    plt.tight_layout()
    print('Done')


def solution_through_time(model, T, node, t_int, dt, time_scale='day', ax=None,
                          **plot_kwargs):
    """Plott solution at a specific node through time
    
    """
    print('Plotting solution through time at node {} '.format(node), end='')
    if ax is None:
        fig, ax = plt.subplots()

    if time_scale == 'day':
        time_factor = 60*60*24
    elif time_scale == 'hour':
        time_factor = 60*60
    else:
        print('Time scale should be hour or day!')

    # t_int in seconds
    number_steps = int(t_int/dt)
    t = np.linspace(0, t_int, number_steps+1)

    ax.plot(t/time_factor, T[node, :], **plot_kwargs)
    ax.set_xlabel('time ({})'.format(time_scale))
    ax.set_ylabel(r'Temperature $^{\circ}C$')
    plt.legend()
    plt.tight_layout()
    print('Done')
    
