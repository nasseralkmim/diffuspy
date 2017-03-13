from diffuspy import boundary
import numpy as np
from diffuspy import backwardeuler
from diffuspy import picard


def solver(model, material, t_int, dt, T0=0, σ_q=None, q_bc=None,
           T_bc=None, T_a=None, h=None, tol=1e-8, plot_reaction_degree=False):
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
    T_np = T0                   # initial Temperature
    α_np = np.zeros(model.ne)   # initial reaction degree
    
    α_avg = [0]
    t_avg = [0]
    # number of defaults parameters in the function
    if len(σ_q.__defaults__) >= 2:
        print('Internal heat source nonlinear on T!')

    # n is the time step, +1 so the loop goes to the last step
    for n in range(1, N+1):
        t = n*dt

        # check nonlinearity
        if len(σ_q.__defaults__) >= 2:
            print('Initializing Iteration at time {:.2f}h with tol: {}!'.format(t/3600, tol))

            # Iteration starts with previous time step converged solution
            # for i=0, T_ip = T_np
            T_ip = T_np

            # iteration for nonlinear case
            i = 0
            error = 1.0
            while error > tol:
                K_u, F_u, α_np = picard.scheme(model, material, dt, t, T_np,
                                               σ_q, q_bc, T_a, h, T_ip, α_np, tol=tol)

                Km, Pm = boundary.temperature(K_u, F_u, model, T_bc)

                T_iu = np.linalg.solve(Km, Pm)

                error = (np.linalg.norm(T_iu - T_ip))/np.linalg.norm(T_iu)

                # update previous iteration solution
                T_ip = T_iu
                
                print('Iteration {} with error {} '.format(i, error))
                i += 1

            T[:, n] = T_iu      # save the converged solution
            print('Solution at time {:.2f}h converged with {} iterations'
                  .format(t/3600, i))

            t_avg.append(t/3600)
            α_avg.append(np.average(α_np))

        else:
            # condensed matrix to form the system K_u T_n+1 = F_u
            K_u, F_u = backwardeuler.scheme(model, material, dt, t, T_np,
                                            σ_q, q_bc, T_a, h)

            Km, Pm = boundary.temperature(K_u, F_u, model, T_bc)

            T[:, n] = np.linalg.solve(Km, Pm)

        # update T previous in time n
        T_np = T[:, n]

    print(' Solution completed!')
    print('Number of elements: {} and Number of nodes: {}'.format(model.ne,
                                                                  model.nn))

    if plot_reaction_degree is True:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot(t_avg, α_avg, 'ks')
        ax.set_xlabel('Time, h')
        ax.set_ylabel(r'Reaction degree $\alpha$')

    print('Temperature field is an array with shape: ', np.shape(T))
    return T
