import numpy as np
from diffuspy.constructor import constructor


def scheme(model, material, dt, t, T_np, σ_q, q_bc, T_a, h, T_ip, α_np, tol=1e-8):
    """perform the backward euler scheme and returns the global K_u and F_u

    The sufix u means update, p means previous

    Args:
        T_np: Previous temperature vector in time n used in backward euler scheme
        T_ip: Previous temperature in the iteration i used to compute α
        α_np: reaction degree from previous time step n used to compute α_up
            using Newton-Raphson scheme

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
        print('{:.2f}h...Building element {}...'.format(t/3600, eid))
        element = constructor(eid, model, material)

        k_q = element.heat_stiffness_matrix(t)
        k_c = element.heat_convection_matrix(h, t)
        k_s = element.heat_capacitance_matrix(t)
        p_t = element.heat_boundary_flux_vector(q_bc, t)
        p_c = element.heat_boundary_convection_vector(T_a, h, t)

        # check if σ_q is non linear dependent on reaction degree
        if 'Reaction Degree' in σ_q.__defaults__:
            α_nu = reaction_degree_iteration(element.conn, T_ip, α_np[eid], dt, t, tol=tol)
            dα = (α_nu - α_np[eid])/dt

            p_q = element.heat_source_vector(σ_q, t, dα=dα)

            # update previous reaction degree in time
            α_np[eid] = α_nu

        elif 'Temperature' in σ_q.__defaults__:
            p_q = element.heat_source_vector(σ_q, t, T_ip)

        K_q[element.id_m] += k_q
        K_c[element.id_m] += k_c
        K_s[element.id_m] += k_s
        P_q[element.id_v] += p_q
        P_t[element.id_v] += p_t
        P_c[element.id_v] += p_c

    F_u = P_q + P_c - P_t + (1/dt)*(K_s @ T_np)
    K_u = (1/dt)*K_s + K_q + K_c

    return K_u, F_u, α_np


def reaction_degree_iteration(conn, T_ip, α_np, dt, t, tol=1e-8):
    """Compute iteratively the reaction degree

    Using the Newton-Raphson scheme

    Args:
        α_np: reaction degree from previous time step n

    Returns:
        α_nu: converged reaction degree

    """
    print('{:.2f}h... Initialize iteration for reaction degree...'.format(t/3600))
    T_avg = np.average(T_ip[conn])

    # data from Kuriakose2016
    nb = 8
    k_n0 = .35e8/3600    # s-1
    A0_k = 1e-5
    Ea_R = 50000/8.314   # K
    α_inf = 0.84

    # number of intermediary steps for the reaction degree
    n = 1000
    dt = dt/n

    # converged solution from previous step for this element
    α_jp = α_np

    for j in range(1, n+1):
        # Iteration in i starts with converged value from previous intermediary step
        α_ip = α_jp

        error = 1.0
        i = 0
        while error > tol:
            df = (- nb/α_inf * k_n0 * (A0_k/α_inf + α_ip) *
                (α_inf - α_ip) * np.exp(- nb*α_ip/α_inf)
                - k_n0 * (A0_k/α_inf + α_ip) * np.exp(- nb*α_ip/α_inf)
                + k_n0 * (α_inf - α_ip) * np.exp(- nb*α_ip/α_inf))
            f = k_n0*(A0_k/α_inf + α_ip)*(α_inf - α_ip)*np.exp(-nb*α_ip/α_inf)

            # convert celsius to kelvin K = C + 273.15
            # Ea_R is in K
            g = np.exp(- Ea_R/(T_avg + 273.15))

            ᴪ = α_ip - dt * f * g - α_jp
            dᴪ = 1 - dt * g * df

            α_iu = α_ip - ᴪ/dᴪ

            error = np.abs(α_iu - α_ip) / np.abs(α_iu)

            # update previous iteration
            α_ip = α_iu
            i += 1

        print('{:.2f}h... Result {} converged with {} iterations'
              .format(t/3600, α_iu, i))

        # converged solution
        α_jp = α_iu

    # after loop in the intermediary steps
    α_ju = α_iu
    
    return α_ju
