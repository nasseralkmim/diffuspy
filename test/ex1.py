from diffuspy import gmsh
from diffuspy import plotter
from diffuspy import data
from diffuspy import transient


model_name = 'patch'

model = gmsh.Parse(model_name)

material = data.Collect()
s = list(model.surf.keys())

material.cndtvt[s[0]] = 1
material.spcfht[s[0]] = 1
material.dnsty[s[0]] = 0.1


def internal_heat(x1, x2, t=1):
    return 10 * x1 * x2 * (10*t)**2


def flux_bc(x1, x2, t=1):
    return {1: 10.0, 2: 0.0, 0: 0.0}


def temperature_bc(x1, x2, t=1):
    return {3: 20.0}

T0 = 10.0

interval = .01
dt = 0.001

frames = transient.solver(model, material, internal_heat, flux_bc,
                          temperature_bc, T0, interval, dt, vmin=0., vmax=100.)


ani = plotter.anime(frames, dt)
ani.save('diffuspy2.gif', writer='imagemagick', bitrate=1000)
# plotter.model(model, ele=True)
plotter.show()
