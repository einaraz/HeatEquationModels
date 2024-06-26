import numpy as np
import matplotlib.pyplot as plt

# This function solves the heat equation in 1D It takes as inputs the length of the rod L, the final time T,
#    the thermal diffusivity alpha, the initial temperature distribution f, the number of spatial grid points Nx, and the number of time steps Nt
#    It returns the temperature distribution at the final time


def heatequation(L, T, alpha, f, Nx, Nt):
    # Initialize the grid
    x   = np.linspace(0, L, Nx+1)    # mesh points in space
    dx  = x[1] - x[0]
    t   = np.linspace(0, T, Nt+1)    # mesh points in time
    dt  = t[1] - t[0]
    F   = alpha*dt/dx**2
    u   = np.zeros(Nx+1)           # unknown u at new time level
    u_n = np.zeros(Nx+1)           # u at the previous time level

    # Set initial condition
    for i in range(0, Nx+1):
        u_n[i] = f(x[i])
    for n in range(0, Nt):
        # Compute u at inner mesh points
        for i in range(1, Nx):
            u[i] = u_n[i] + F*(u_n[i-1] - 2*u_n[i] + u_n[i+1])
        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        # Update u_n before next step
        u_n[:]= u

    return u

# Define the initial temperature distribution
def f(x):
    return np.sin(np.pi*x)

if __name__ == '__main__':
    # Plot the solution for a specific set of parameters
    u = heatequation(L=1, T=0.1, alpha=0.1, f=f, Nx=50, Nt=50)
    plt.plot(np.linspace(0, 1, 51), u)
    plt.show()
