import numpy as np
from diffuspy.constructor import constructor


def scheme(model, material, dt, t, T_p, σ_q, q_bc, T_a):
    """perform the backward euler scheme and returns the global K_u and F_u

    The sufix u means update, p means previous

    Return:
        K_u: the equivalent stiffness matrix updated
        F_u: the equivalent load vector updated

    """
    K_q = np.zeros((model.ndof, model.ndof))
    K_c = np.zeros((model.ndof, model.ndof))
    K_s = np.zeros((model.ndof, model.ndof))
    P_q = np.zeros(model.ndof)
    P_t = np.zeros(model.ndof)
    P_c = np.zeros(model.ndof)

    for eid, type in enumerate(model.TYPE):
        element = constructor(eid, model, material)

        k_q = element.heat_stiffness_matrix(t)
        k_c = element.heat_convection_matrix(t)
        k_s = element.heat_capacitance_matrix(t)
        p_q = element.heat_source_vector(σ_q, t)
        p_t = element.heat_boundary_flux_vector(q_bc, t)
        p_c = element.heat_boundary_convection_vector(T_a, t)

        K_q[element.id_m] += k_q
        K_c[element.id_m] += k_c
        K_s[element.id_m] += k_s
        P_q[element.id_v] += p_q
        P_t[element.id_v] += p_t
        P_c[element.id_v] += p_c

    F_u = P_q + P_c - P_t + (1/dt)*(K_s @ T_p)
    K_u = (1/dt)*K_s + K_q + K_c

    return K_u, F_u
