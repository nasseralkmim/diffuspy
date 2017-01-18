from diffuspy.mesh import gmsh
from diffuspy.postprocess import plotter
from diffuspy.material import Material
from diffuspy.solvers import steadystate

model_name = 'patch'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(λ={s[0]: 1},
                    c={s[0]: 1},
                    ρ={s[0]: 0.1})

def σ_q(x1, x2, t=1):
    return 0

def q_bc(x1, x2, t=1):
    return {1: -20}

def T_bc(x1, x2, t=1):
    return {3: 50}

T = steadystate.solver(model, material, σ_q , q_bc, T_bc)

plotter.contour(model, T)
plotter.model(model, ele=True, edges_label=True)

plotter.show()
