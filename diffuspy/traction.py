import numpy as np


def Pt_vector(model, f, t=1):
    """Build the global vector due boundary traction

    """
    Pt = np.zeros(model.nn)

    for bound_line in f(1, 1).keys():
        for e, side, line in model.bound_ele:
            if line == bound_line:
                conn = model.CONN[e]
                xyz = model.XYZ[conn]
                pt = pt_vector(model, xyz, side, line, f, t)

                id = conn

                Pt[id] += pt

    return Pt


def pt_vector(model, xyz, side, line, f, t=1):
    """Build element load vector due traction boundary condition

    """
    pt = np.zeros(4)

    gp = np.array([
        [[-1.0/np.sqrt(3), -1.0],
         [1.0/np.sqrt(3), -1.0]],
        [[1.0, -1.0/np.sqrt(3)],
         [1.0, 1.0/np.sqrt(3)]],
        [[-1.0/np.sqrt(3), 1.0],
         [1.0/np.sqrt(3), 1.0]],
        [[-1.0, -1.0/np.sqrt(3)],
         [-1.0, 1/np.sqrt(3)]]])

    for w in range(2):
        model.basis_function(gp[side, w])
        model.jacobian(xyz)
        dL = model.ArchLength[side]
        x1, x2 = model.mapping(xyz)

        pt += model.phi[:]*f(x1, x2, t)[line]*dL

    return pt
