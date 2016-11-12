from diffuspy import boundconditions
from diffuspy import stiffness
from diffuspy import load
from diffuspy import traction
import numpy as np


def solver(model, material, b_heat, flux_bc,
           temp_bc, t=1):

    K = stiffness.K_matrix(model, material)

    Pq = load.Pq_vector(model, b_heat, t)

    Pt = traction.Pt_vector(model, flux_bc, t)

    P = Pq + Pt

    Kf, Pf = boundconditions.temperature(K, P, model, temp_bc, t)

    T = np.linalg.solve(Kf, Pf)

    return T
