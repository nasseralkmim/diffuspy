from pyisson import fem
from pyisson import gmsh
from pyisson import plotter
from pyisson import data

meshName = 'square'

mesh = gmsh.parse(meshName)

material = {'k-c-rho': [10.0, 1.0, 0.1]}
surf = data.Collect()
surf.k[0] = 10
surf.c[0] = 1
surf.rho[0] = 0.1


def internal_heat(x1, x2, t=1):
    return x1*x2*t


def flux_bc(x1, x2):
    return {
        1: 100.0,
        2: 0.0,
        0: 0.0
    }


def temperature_bc(x1, x2):
    return {
        3: 10.0,
    }

T = fem.solver(mesh, material, internal_heat, flux_imposed,
               temperature_imposed)

plotter.contour(mesh, T, cmap='hot', lev=10)
plotter.model(mesh, ele=True)

plotter.show()
