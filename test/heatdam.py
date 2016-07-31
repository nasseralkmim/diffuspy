from diffuspy import gmsh
from diffuspy import plotter
from diffuspy import data
from diffuspy import transient
from diffuspy import steadystate
import time


start_time = time.time()

model_name = 'barragem3'

model = gmsh.Parse(model_name)

material = data.Collect()
s = list(model.surf.keys())

material.cndtvt[s[0]] = 100
material.spcfht[s[0]] = 0.05
material.dnsty[s[0]] = 10

material.cndtvt[s[1]] = 1
material.spcfht[s[1]] = 0.03
material.dnsty[s[1]] = 10


def internal_heat(x1, x2, t=1):
    return (100*t)**2


def flux_bc(x1, x2, t=1):
    return {7: 100, 8: 100, 9: 100, 10: 100}


def temperature_bc(x1, x2, t=1):
    return {2: 20, 3: 20, 4: 20, 5: 20, 6: 20, 0: 5, 1: 5}

T0 = 5.0

interval = 0.1
dt = 0.001

frames = transient.solver(model, material, internal_heat, flux_bc,
                          temperature_bc, T0, interval, dt, vmin=0, vmax=40)

ani = plotter.anime(frames, dt)
ani.save('2surf_dam8.gif', writer='imagemagick', bitrate=1000)

T = steadystate.solver(model, material, internal_heat, flux_bc,
                       temperature_bc)

# plotter.contour(model, T)
# plotter.model(model, edges_label=True, surf=True, surf_label=True)
# plotter.model(model, ele=True)
plotter.show()

print(interval/dt * (model.nn**2 + model.nn**3) * 100 / (1.9*1e9))

print("--- %s seconds ---" % (time.time() - start_time))
