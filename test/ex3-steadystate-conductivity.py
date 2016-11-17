from diffuspy import steadystate
from diffuspy import gmsh
from diffuspy import plotter
from diffuspy.material import Material

model_name = 'ex3-steady'

model = gmsh.Parse(model_name)

h1, h2 = .3, .5
d1, d2 = .4, 1
k1, k2 = 100, 5

def conductivity(x1, x2, t=1):
    return ((h1 <= x1 <= h1+d1 and h2 <= x2 <= h2+d2)*k1 +
            (x1 < h1 or x2 < .5)*k2 +
            (x1 > h1+d1 or x2 > h2+d2)*k2)

s = list(model.surf.keys())
material = Material(cndtvt={s[0]: conductivity},
                    spcfht={s[0]: 1},
                    dnsty={s[0]: 0.1})

def internal_heat(x1, x2, t=1):
    return 0.0

def flux_bc(x1, x2, t=1):
    return {}

def temperature_bc(x1, x2, t=1):
    return {1: 32, 3: 50}

T = steadystate.solver(model, material, internal_heat, flux_bc,
                       temperature_bc)

plotter.contour(model, T)
# import matplotlib.pyplot as plt
# plt.savefig('ex-heat-st-1.pdf', bbox_inches='tight')
plotter.model(model, ele=True, edges_label=True)
# plt.savefig('ex-heat-st-2.pdf', bbox_inches='tight')
plotter.show()
