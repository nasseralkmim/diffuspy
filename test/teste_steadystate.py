from diffuspy.mesh import gmsh
from diffuspy.postprocess import plotter
from diffuspy.material import Material
from diffuspy.solvers import steadystate

model_name = 'patch-refined'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(λ={s[0]: 1},  # condutivity W/m C
                    c={s[0]: 1},  # specific heat J/kg C
                    ρ={s[0]: 0.1})  # density kg/ m3


# Internal heat (q) current density source (σ) in W/ m3
def σ_q(x1, x2, t=1, T=1):
    return 0


# Heat (q) current density assigned at the boundary (bc) in W/m3
def q_bc(x1, x2, t=1):
    return {}


# Temperature (T) assigned at boundary (bc) in K or C
def T_bc(x1, x2, t=1):
    return {3: 30}


# Surface condutance in W/m2 C
def h(x1, x2, t=1):
    return {1: 100}


# Temperature of the air K or C
def T_a(x1, x2, t=1):
    return {1: 1000}


T = steadystate.solver(model, material,
                       σ_q=σ_q, q_bc=q_bc, T_bc=T_bc, T_a=T_a, h=h)

plotter.contour(model, T)
plotter.model(model, ele=True, edges_label=True, nodes_label=True)

plotter.show()
