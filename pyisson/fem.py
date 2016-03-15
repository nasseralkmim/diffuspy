from scipy.sparse.linalg import spsolve
from scipy import sparse
from pyisson import element1dof
from pyisson import assemble1dof
from pyisson import boundaryconditions1dof


def solver(mesh, material, internal_heat, flux_imposed,
           temperature_imposed):

    ele = element1dof.Matrices(mesh)

    # Surfaces tags
    s = mesh.surfaces

    if len(s) < len(material):
        raise ValueError(
            'There are more Materials assigned than Surfaces defined!')

    if len(s) > len(material):
        raise ValueError(
            'There are more Surfaces defined than Material assigned!')

    mat = {s[i]: material[j] for i, j in enumerate(material)}

    ele.stiffness(mat)

    ele.load(internal_heat)

    K = assemble1dof.globalMatrix(ele.K, mesh)

    P0q = assemble1dof.globalVector(ele.R, mesh)

    P0t = boundaryconditions1dof.neumann(mesh, flux_imposed)

    P0 = P0q + P0t

    K, P0 = boundaryconditions1dof.dirichlet(K, P0, mesh, temperature_imposed)

    Ks = sparse.csc_matrix(K)

    U = spsolve(Ks, P0)

    return U
