import numpy as np


def Pq_vector(model, q, t=1):
    """Build the load vector for the internal heat source

    """
    Pq = np.zeros(model.nn)
    for e, conn in enumerate(model.CONN):
        xyz = model.XYZ[conn]
        pq = pq_vector(model, xyz, q, t)

        id = conn

        Pq[id] += pq

    return Pq


def pq_vector(model, xyz, q, t=1):
    """Build the element vector for the internal heat source

    """
    gauss_points = model.chi / np.sqrt(3.0)

    pq = np.zeros(4)
    for gp in gauss_points:
        model.basis_function(gp)
        model.jacobian(xyz)
        dJ = model.detJac
        x1, x2 = model.mapping(xyz)

        pq[:] += model.phi[:]*q(x1, x2, t)*dJ

    return pq
