from scipy.sparse.linalg import spsolve
from scipy import sparse
from diffuspy import boundconditions
from diffuspy import stiffness
from diffuspy import load
from diffuspy import traction


def solver(model, material, internal_heat, flux_bc,
           temperature_bc, t=1):

    K = stiffness.K_matrix(model, material)

    Pq = load.Pq_vector(model, internal_heat, t)

    Pt = traction.Pt_vector(model, flux_bc, t)

    P = Pq + Pt

    Kf, Pf = boundconditions.temperature(K, P, model, temperature_bc, t)

    Ks = sparse.csc_matrix(K)

    T = spsolve(Ks, Pf)

    return T
