__author__ = 'Nasser'

import heat2dtransient

meshName = 'square'

material= {'k-c-rho':[1.0, 1.0, 0.1]}

def internal_heat(x1, x2, t):
    return 0.0

def flux_imposed(x1, x2, t):
    return {
        1:10.0,
        2:0.0,
        0:0.0
    }

def temperature_imposed(x1, x2, t):
    return {
            3:20.0,
            }

delta_t = 0.0001
time_interval = 1001

heat2dtransient.solver(meshName,
                       delta_t,
                       time_interval,
                       material,
                       internal_heat,
                       flux_imposed,
                       temperature_imposed,
                       plot_at=[0.0, 0.001, 0.00001, 0.01, 0.0001])

