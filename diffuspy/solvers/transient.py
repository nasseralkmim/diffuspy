from diffuspy import boundary
import numpy as np
from diffuspy import backwardeuler


def solver(model, material, t_int, dt, T0=0, σ_q=None, q_bc=None,
           T_bc=None, T_a=None, h=None):
    """Solver for the transient problem

    """
    print('Initializing solver...')
    # compute the number of time steps N
    # dt=t_int/3 -> N=3, t_int=6 -> 0-2, 2-4, 4-6
    N = int(t_int/dt)

    print('Number of steps: ', N)
    # plus 1 to consider T[t=0]
    T = np.zeros((model.ndof, N+1))

    # initial condition
    if np.size(T0) == 1:
        T0 = np.ones(model.ndof) * T0

    T[:, 0] = T0
    T_p = T0

    # n is the time step, +1 so the loop goes to the last step
    for n in range(1, N+1):
        t = n*dt

        K_u, F_u = backwardeuler.scheme(model, material, dt, t, T_p,
                                        σ_q, q_bc, T_a, h)

        Km, Pm = boundary.temperature(K_u, F_u, model, T_bc)

        T[:, n] = np.linalg.solve(Km, Pm)

        T_p = T[:, n]           # update T previous

    print('Solution completed!')
    print('Temperature field is an array with shape: ', np.shape(T))
    return T
