import numpy as np


def M_matrix(model, material):
    """Build the global mass matrix

    """
    M = np.zeros((model.nn, model.nn))
    for e, conn in enumerate(model.CONN):
        xyz = model.XYZ[conn]
        surf = model.surf_of_ele[e]
        dnsty = material.dnsty[surf]
        spcfht = material.spcfht[surf]

        m = m_matrix(model, xyz, dnsty, spcfht)
        id = np.ix_(conn, conn)

        M[id] += m

    return M


def m_matrix(model, xyz, dnsty, spcfht):
    """Build element mass matrix

    """
    gauss_points = model.chi / np.sqrt(3.0)
    m = np.zeros((4, 4))
    for gp in gauss_points:
        model.basis_function(gp)
        model.jacobian(xyz)
        dJ = model.detJac

        N = np.array([model.phi])

        m += N.T @ N * dJ * dnsty * spcfht

    return m
