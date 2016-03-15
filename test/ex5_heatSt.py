from pyisson import fem
from pyisson import gmsh
from pyisson import plotter

meshName = 'square'

mesh = gmsh.parse(meshName)

material = {'k-c-rho': [10.0, 1.0, 0.1]}


def internal_heat(x1, x2):
    return 0.0


def flux_imposed(x1, x2):
    return {
        1: 100.0,
        2: 0.0,
        0: 0.0
    }


def temperature_imposed(x1, x2):
    return {
        3: 10.0,
    }

T = fem.solver(mesh, material, internal_heat, flux_imposed,
               temperature_imposed)

plotter.contour(mesh, T, cmap='hot', lev=10)
plotter.model(mesh, ele=True)

plotter.show()
