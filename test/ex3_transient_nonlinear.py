from diffuspy.mesh import gmsh
from diffuspy.postprocess import plotter
from diffuspy.material import Material
from diffuspy.solvers import transient
import numpy as np

model_name = 'patch-refined'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(λ={s[0]: 1},  # condutivity W/m C
                    c={s[0]: 1},  # specific heat J/kg C
                    ρ={s[0]: 1})  # density kg/ m3


# Internal heat (q) current density source (σ) in W/ m3
def σ_q(x1, x2, t=1, T=1):
    return T**2*(t/3600)*1e-3


# Heat (q) current density assigned at the boundary (bc) in W/m3
# negative mean entrying the domain
def q_bc(x1, x2, t=1):
    return {}


# Temperature (T) assigned at boundary (bc) in K or C
def T_bc(x1, x2, t=1):
    return {3: 30, 1: 30}


t_int = 60 * 60 * 10
dt = t_int / 100
T0 = 30

T = transient.solver(model, material, t_int, dt,
                     T0, σ_q, q_bc, T_bc)

plotter.model(model, ele=True, nodes_label=True)
plotter.array(T, node=36, interval=t_int, dt=dt)
plotter.contour_animation(model, T, t_int, dt, interval=200,
                          time_text_color='black', name='temp_field.gif')
plotter.show()
