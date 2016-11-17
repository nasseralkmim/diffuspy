from diffuspy import steadystate
from diffuspy import gmsh
from diffuspy import plotter
from diffuspy.material import Material

model_name = 'patch'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(cndtvt={s[0]: 1},
                    spcfht={s[0]: 1},
                    dnsty={s[0]: 0.1})

def internal_heat(x1, x2, t=1):
    return 0.0


def flux_bc(x1, x2, t=1):
    return {1: 10.0, 2: 0.0, 0: 0.0}


def temperature_bc(x1, x2, t=1):
    return {3: 10.0}

T = steadystate.solver(model, material, internal_heat, flux_bc,
                       temperature_bc)

plotter.contour(model, T)
plotter.model(model, ele=True, edges_label=True)

plotter.show()
