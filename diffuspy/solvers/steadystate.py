from diffuspy import boundary
import numpy as np
from diffuspy.constructor import constructor


def solver(model, material, σ_q=None, q_bc=None,
           T_bc=None, T_a=None, h=None):
    """Solver for the steadystate problem

    """
    print('Initializing solver...', end='')
    K_q = np.zeros((model.ndof, model.ndof))
    K_c = np.zeros((model.ndof, model.ndof))
    P_q = np.zeros(model.ndof)
    P_t = np.zeros(model.ndof)
    P_c = np.zeros(model.ndof)

    for eid, type in enumerate(model.TYPE):
        element = constructor(eid, model, material)

        k_q = element.heat_stiffness_matrix()
        k_c = element.heat_convection_matrix(h)
        p_q = element.heat_source_vector(σ_q)
        p_t = element.heat_boundary_flux_vector(q_bc)
        p_c = element.heat_boundary_convection_vector(T_a, h)

        K_q[element.id_m] += k_q
        K_c[element.id_m] += k_c
        P_q[element.id_v] += p_q
        P_t[element.id_v] += p_t
        P_c[element.id_v] += p_c

    P = P_q + P_c - P_t
    K = K_q + K_c

    Km, Pm = boundary.temperature(K, P, model, T_bc)

    T = np.linalg.solve(Km, Pm)
    print('Solution completed!')
    return T
