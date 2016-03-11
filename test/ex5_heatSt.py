__author__ = 'Nasser'

import poisson.heat2d as heat2d

meshName = 'square'

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

heat2d.solver(meshName,
              material,
              internal_heat,
              flux_imposed,
              temperature_imposed,
              plotUndeformed={'Domain': True,
                              'Elements': True,
                              'NodeLabel': False,
                              'EdgesLabel':  True,
                              'ElementLabel': False,
                              'SurfaceLabel': False},
              plotTemperature={'Contour': True}
              )


