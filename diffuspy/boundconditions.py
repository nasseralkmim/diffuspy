def temperature(K, P, model, T_bc, t=1):
    """Apply Dirichlet boundary conditions

    """
    for bound_line in T_bc(1, 1).keys():
        for line, n1, n2 in model.nodes_in_bound_line:
            if line == bound_line:
                K[n1, :] = 0.0
                K[n2, :] = 0.0

                K[n1, n1] = 1.0
                K[n2, n2] = 1.0

                t1 = T_bc(model.XYZ[n1, 0], model.XYZ[n1, 1], t)
                t2 = T_bc(model.XYZ[n2, 0], model.XYZ[n2, 1], t)

                P[n1] = t1[line]
                P[n2] = t2[line]

    return K, P
