from diffuspy import boundconditions
from diffuspy import stiffness
from diffuspy import load
from diffuspy import traction
import numpy as np
from diffuspy.constructor import constructor

def solver(model, material, body_heat, flux_bc,
           temperature_bc, t=1):
    """Solver for the steadystate problem

    """
    Kq = np.zeros((model.ndof, model.ndof))
    Pq = np.zeros(model.ndof)
    Pt = np.zeros(model.ndof)

    for eid, type in enumerate(model.TYPE):
        element = constructor(eid, model, material)

        k = element.stiffness_matrix()
        pq = element.heat_body_vector(body_heat)
        pt = element.heat_bound_vector(flux_bc)

        Kq[element.id_m] += k
        Pq[element.id_v] += pq
        Pt[element.id_v] += pt

    P = Pq + Pt

    Km, Pm = boundary.dirichlet(Kq, P, model, temperature_bc)

    T = np.linalg.solve(Km, Pm)

    return T
