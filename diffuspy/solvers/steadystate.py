from diffuspy import boundary
from diffuspy import stiffness
from diffuspy import load
from diffuspy import traction
import numpy as np
from diffuspy.constructor import constructor

def solver(model, material, σ_q=None, q_bc=None,
           T_bc=None, t=1):
    """Solver for the steadystate problem

    """
    Kq = np.zeros((model.ndof, model.ndof))
    Pq = np.zeros(model.ndof)
    Pt = np.zeros(model.ndof)

    for eid, type in enumerate(model.TYPE):
        element = constructor(eid, model, material)

        k = element.heat_stiffness_matrix()
        pq = element.heat_source_vector(σ_q)
        pt = element.heat_boundary_vector(q_bc)

        Kq[element.id_m] += k
        Pq[element.id_v] += pq
        Pt[element.id_v] += pt

    P = Pq + Pt

    Km, Pm = boundary.temperature(Kq, P, model, T_bc)

    T = np.linalg.solve(Km, Pm)

    return T
