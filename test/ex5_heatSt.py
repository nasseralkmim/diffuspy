from diffuspy.mesh import gmsh
from diffuspy.postprocess import plotter
from diffuspy.material import Material
from diffuspy.solvers import steadystate

model_name = 'patch'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(λ={s[0]: 1}, # condutivity W/m K
                    c={s[0]: 1}, # specific heat J/kg K
                    ρ={s[0]: 0.1}) # density kg/ m3

# Internal heat (q) current density source (σ) in W/ m3
def σ_q(x1, x2, t=1):
    return 0

# Heat (q) current density assigned at the boundary (bc) in W/m3
def q_bc(x1, x2, t=1):
    return {1: -20}

# Temperature (T) assigned at boundary (bc) in K or C
def T_bc(x1, x2, t=1):
    return {3: 50}

T = steadystate.solver(model, material, σ_q , q_bc, T_bc)

plotter.contour(model, T)
plotter.model(model, ele=True, edges_label=True)

plotter.show()
