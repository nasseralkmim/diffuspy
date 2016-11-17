import numpy as np


def K_matrix(model, material, t=1):
    """Build the global stiffness matrix

    """
    K = np.zeros((model.nn, model.nn))
    for e, conn in enumerate(model.CONN):

        xyz = model.XYZ[conn]
        surf = model.surf_of_ele[e]

        try:
            cndtvt = material.cndtvt[surf]
        except:
            print('Surface {} has no property assigned!'
                  'Default 1.0 was used!'.format(surf))
            cndtvt = 1.0

        k = k_matrix(model, xyz, cndtvt, t)

        id = np.ix_(conn, conn)

        K[id] += k

    return K


def k_matrix(model, xyz, cndtvt, t):
    """Build the element stiffness matrix

    """
    gauss_points = model.chi / np.sqrt(3.0)

    k = np.zeros((4, 4))
    for gp in gauss_points:
        model.basis_function(gp)
        model.jacobian(xyz)
        dJ = model.detJac

        B = model.dphi_xi

        # Check if condutivity is a functions 
        if callable(cndtvt) is True:
            x1, x2 = model.mapping(xyz)
            cndtvt_value = cndtvt(x1, x2, t)
        else:
            cndtvt_value = cndtvt

        k += cndtvt_value*(B.T @ B)*dJ

    return k
