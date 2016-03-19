from diffuspy import steadystate
from diffuspy import gmsh
from diffuspy import plotter
from diffuspy import data

model_name = 'patch'

model = gmsh.Parse(model_name)

material = data.Collect()
s = list(model.surf.keys())

material.cndtvt[s[0]] = 1
material.spcfht[s[0]] = 1
material.dnsty[s[0]] = 0.1


def internal_heat(x1, x2, t=1):
    return 0.0


def flux_bc(x1, x2, t=1):
    return {1: 10, 2: 0.0, 0: 0.0}


def temperature_bc(x1, x2, t=1):
    return {3: 10.0}

T = steadystate.solver(model, material, internal_heat, flux_bc,
                       temperature_bc)

plotter.contour(model, T)
plotter.model(model, ele=True, edges_label=True)

plotter.show()
