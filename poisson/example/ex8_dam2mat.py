__author__ = 'Nasser'

import numpy as np
import elasticity2d

meshName = 'barragem3'

material = {'Ef-nu': [30e9, 0.17], 'Eb-nu': [15e9, 0.17]}

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
                                 'strain':False}
                    )


