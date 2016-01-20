__author__ = 'Nasser'

import heat2d

meshName = 'barragem3'

material= {'k-c-rho-1':[2.0, 800.0, 2000.0], 'k-c-rho-2':[1.3, 1.0, 2000.0]}

def internal_heat(x1, x2):
    return (x2 > 70.0)*0.4

def flux_imposed(x1, x2):
    return {
    }

def temperature_imposed(x1, x2):
    return {
            0:10.0,
            1:10.0,
            2:30.0,
            3:30.0,
            4:30.0,
            5:30.0,
            6:30.0,
            7:30.0,
            8:15.0,
            9:15.0,
            10:10.0
            }

heat2d.solver(meshName,
              material,
              internal_heat,
              flux_imposed,
              temperature_imposed,
              plotUndeformed={'Domain':True,
                              'Elements':True,
                              'NodeLabel':False,
                              'EdgesLabel':False,
                              'ElementLabel':False,
                              'SurfaceLabel':False},
              plotTemperature={'Contour':True}
              )


