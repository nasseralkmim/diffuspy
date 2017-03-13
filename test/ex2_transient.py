from diffuspy.mesh import gmsh
from diffuspy.postprocess import plotter
from diffuspy.material import Material
from diffuspy.solvers import transient
import numpy as np

model_name = 'patch'

model = gmsh.Parse(model_name)

s = list(model.surf.keys())
material = Material(λ={s[0]: 1},  # condutivity W/m C
                    c={s[0]: 1},  # specific heat J/kg C
                    ρ={s[0]: 1})  # density kg/ m3


# Internal heat (q) current density source (σ) in W/ m3
def σ_q(x1, x2, t=1):
    return 0


# Heat (q) current density assigned at the boundary (bc) in W/m3
# negative mean entrying the domain
def q_bc(x1, x2, t=1):
    return {}


# Temperature (T) assigned at boundary (bc) in K or C
def T_bc(x1, x2, t=1):
    return {3: 30}



# Surface condutance in W/m2 C
def h(x1, x2, t=1):
    return {1: 100}


# Temperature of the air in C
def T_a(x1, x2, t=1):
    return {1: 100 * np.sin(2*np.pi*t/(60*60))}


t_int = 60 * 60 * 10
dt = t_int / 10
T0 = 30

plotter.function(T_a, t_int, line=1, xlabel='Time, (h)',
                 ylabel=r'Temperature, ($^{\circ}C$)', c='black', marker='x')

T = transient.solver(model, material, t_int, dt,
                     T0, σ_q, q_bc, T_bc, T_a, h)

plotter.contour_animation(model, T, t_int, dt, interval=200,
                          time_text_color='black', name='temp_fiel.gif')
plotter.show()
