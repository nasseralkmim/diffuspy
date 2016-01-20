__author__ = 'Nasser'

import numpy as np
import elasticity2d
import heat2d

meshName = 'barragem3'

materialHeat= {'k-c-rho-1':[2.0, 800.0, 2000.0], 'k-c-rho-2':[1.3, 1.0, 2000.0]}

def internal_heat(x1, x2):
    return 0.0

def flux_imposed(x1, x2):
    return {
    }

airT = 35.0
waterT = 15.0
foundT = 5.0
def temperature_imposed(x1, x2):
    return {
            0:foundT,
            1:foundT,
            2:airT,
            3:airT,
            4:airT,
            5:airT,
            6:airT,
            7:airT,
            8:waterT,
            9:waterT,
            10:waterT
            }

T2 = heat2d.solver(meshName,
              materialHeat,
              internal_heat,
              flux_imposed,
              temperature_imposed,
              plotUndeformed={'Domain':False,
                              'Elements':False,
                              'NodeLabel':False,
                              'EdgesLabel':False,
                              'ElementLabel':False,
                              'SurfaceLabel':False},
              plotTemperature={'Contour':False}
              )


materialHeat= {'k-c-rho-1':[1.8, 800.0, 2000.0], 'k-c-rho-2':[1.3, 1.0, 2000.0]}

def internal_heat(x1, x2):
    return 0.0

def flux_imposed(x1, x2):
    return {
    }

airT = 20.0
waterT = 12.0
foundT = 5.0
def temperature_imposed(x1, x2):
    return {
            0:foundT,
            1:foundT,
            2:airT,
            3:airT,
            4:airT,
            5:airT,
            6:airT,
            7:airT,
            8:waterT,
            9:waterT,
            10:waterT
            }

T1 = heat2d.solver(meshName,
              materialHeat,
              internal_heat,
              flux_imposed,
              temperature_imposed,
              plotUndeformed={'Domain':False,
                              'Elements':False,
                              'NodeLabel':False,
                              'EdgesLabel':False,
                              'ElementLabel':False,
                              'SurfaceLabel':False},
              plotTemperature={'Contour':False}
              )

material = {'Ef-nu': [16e9, 0.25], 'Eb-nu': [24e9, 0.17]}

dT = T2 - T1

def body_forces(x1, x2): 
    return np.array([
        0.0,
        (x2 >= 70.0)*-24e3+(x2 < 70.0)*-30e3,
    ])

gma = 9820.0
def traction_imposed(x1, x2): 
    return {
        ('line', 7): [(x2 <= 186.13)*(gma*(186.13-x2)), 0.0],
        ('line', 8): [(17.07*gma - 116.13*gma)*(x2-70.0)/(169.06-70.0) +
                                                116.13*gma,
                      -((17.07*gma - 116.13*gma)*(x2-70.0)/(169.06-70.0) +
                        116.13*gma)],
        ('line', 9): [0.0, - 116.13*gma]
    }

def displacement_imposed(x1, x2): 
    return {
        ('line', 0): [0.0, 0.0],
        ('line', 10): [0.0, 'free'],
        ('line', 1): [0.0, 'free']
    }


elasticity2d.solver(meshName,
                    material,
                    body_forces,
                    traction_imposed,
                    displacement_imposed,
                    plotUndeformed={'Domain': True,
                                    'Elements': False,
                                    'NodeLabel': False,
                                    'EdgesLabel': True,
                                    'ElementLabel': False,
                                    'SurfaceLabel': True},
                    plotStress={'s11': False,
                                's22': False,
                                's12': False,
                                'sPmax': True,
                                'sPmin': True},
                    plotDeformed={'DomainUndeformed': True,
                                  'ElementsUndeformed': False,
                                  'DomainDeformed': True,
                                  'ElementsDeformed': True,
                                  'DeformationMagf': 1000},
                    printReport={'U':False,
                                 'stress':False,
                                 'strain':False},
                    initT=dT
                    )


