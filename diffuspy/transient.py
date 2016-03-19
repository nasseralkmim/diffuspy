from diffuspy import boundconditions
from diffuspy import stiffness
from diffuspy import load
from diffuspy import traction
from diffuspy import mass
from diffuspy import plotter
from scipy.linalg import solve
import numpy as np
import matplotlib.pylab as plt


def solver(model, material, internal_heat, flux_bc,
           temperature_bc, T0, interval, dt, vmin=None, vmax=None):

    K = stiffness.K_matrix(model, material)

    M = mass.M_matrix(model, material)

    Tn = np.full(model.nn, T0)

    frames = []

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, .8, .8])
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_title(r'Temperature ($^{\circ}C$)')

    cs = plotter.contour(model, Tn, contour_label=False, cbar=False,
                         figure=False, vmin=vmin, vmax=vmax)

    im = cs.collections

    frames.append(im)
    vn, vx = 0, 0
    for step in range(1, int(interval/dt + 1)):
        t = step * dt

        Pq = load.Pq_vector(model, internal_heat, t)
        Pt = traction.Pt_vector(model, flux_bc, t)
        P = Pq + Pt

        W = dt * (P - K @ Tn) + M @ Tn
        Mf, Wf = boundconditions.temperature(M, W, model, temperature_bc, t)
        Tn = solve(Mf, Wf)

        if np.amin(Tn) < vn:
            vn = np.amin(Tn)

        if np.amax(Tn) > vx:
            vx = np.amax(Tn)

        cs = plotter.contour(model, Tn, contour_label=False,
                             cbar=False, figure=False, vmin=vmin, vmax=vmax)

        im = cs.collections

        frames.append(im)

    print(vn, vx)

    # Change the colorbar range
    sm = plt.cm.ScalarMappable(cmap='hot', norm=plt.Normalize(vmin=vmin,
                                                              vmax=vmax))
    # fake up the array of the scalar mappable. Urgh...
    sm._A = []
    plt.colorbar(sm)
    
    return frames
