from diffuspy import gmsh
from diffuspy import plotter
from diffuspy import data
from diffuspy import transient
from diffuspy import steadystate


model_name = 'mesh2surf2'

model = gmsh.Parse(model_name)

material = data.Collect()
s = list(model.surf.keys())


material.cndtvt[s[0]] = 1
material.spcfht[s[0]] = 1
material.dnsty[s[0]] = 0.1

material.spcfht[s[1]] = 1
material.dnsty[s[1]] = 0.1
material.cndtvt[s[1]] = 10


def internal_heat(x1, x2, t=1):
    return 0.0


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

# T = steadystate.solver(model, material, internal_heat, flux_bc,
#                        temperature_bc)

# plotter.contour(model, T)

# plotter.model(model, edges_label=True, surf=True, surf_label=True)
# plotter.show()
