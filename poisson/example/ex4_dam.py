__author__ = 'Nasser'

import numpy as np
import elasticity2d

meshName = 'barragem2'

material = {'E-nu': [34473.786*1e6, 0.1]}

def body_forces(x1, x2): 
    return np.array([
        0.0,
        -24357.0,
    ])

gma = 9820.0*10
def traction_imposed(x1, x2): 
    return {
        ('line', 7): [(x2 <= 186.13)*(-gma*(x2-186.13)), 0.0],
        ('line', 8): [(17.07*gma - 186.13*gma)*x2/169.06 + 186.13*gma,
                      -((17.07*gma - 186.13*gma)*x2/169.06 + 186.13*gma)],
        ('line', 9): [0.0, - 186.13*gma]
    }

def displacement_imposed(x1, x2): 
    return {
    ('line', 0): [0.0, 0.0],
    ('line', 1): [0.0, 0.0],
    ('line', 10): [0.0, 0.0]
    }

elasticity2d.solver(meshName,
                    material,
                    body_forces,
                    traction_imposed,
                    displacement_imposed,
                    plotUndeformed={'Domain': True,
                                    'Elements': True,
                                    'NodeLabel': False,
                                    'EdgesLabel': False,
                                    'ElementLabel': False,
                                    'SurfaceLabel':False},
                    plotStress={'s11': False,
                                's22': False,
                                's12': False,
                                'sPmax': True,
                                'sPmin': True},
                    plotDeformed={'DomainUndeformed': True,
                                  'ElementsUndeformed': True,
                                  'DomainDeformed': True,
                                  'ElementsDeformed': True,
                                  'DeformationMagf': 500},
                    printReport={'U':False,
                                 'stress':False,
                                 'strain':False}
                    )


