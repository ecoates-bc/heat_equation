# Solve the 1D heat equation with a system of ODEs
# Taken from https://youtu.be/4Bfjs7yemJ8?list=PLcTpn5-ROA4xkuVzXOqIccUVyhzGUN9dw&t=143

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib.animation import FuncAnimation

L = 1  # length of a rod
T = 0.05  # total time
k = 1  # heat constant
N = 50  # samples on the rod
M = 100  # timesteps

dx = L / N
dt = T / M

# create matrix A
A = np.zeros((N-2, N-2))
for n in range(0, N-2):
    A[n, n] = -2 * 1/dx**2
    if n > 0:
        A[n, n-1] = 1/dx**2
    if n < N-3:
        A[n, n+1] = 1/dx**2

# At x=0 and x=N, the initial values are zero, so we can leave out a b vector


# Initial condition
def position_at(i):
    return i * dx


u0 = np.zeros((N-2))
for i in range(N-2):
    u0[i] = math.sin(4 * math.pi * position_at(i))


# define the ODE system
def ode_system_func(u, t):
    du_dt = np.matmul(A, u)
    return du_dt


# define time samples
t = np.linspace(0, T, M)

# solve
soln = odeint(ode_system_func, u0, t)

# plot an animation
fig = plt.figure()
ax = plt.axes(xlim=(0,L), ylim=(-1,1))
plt.style.use("seaborn-pastel")
line, = ax.plot([], [], lw=3)


def animate(i):
    x = np.linspace(0, L, N)
    if i == 0 or i == N-1:
        y = np.zeros(N)
    else:
        y = np.concatenate(([0], soln[i,:], [0]))
    line.set_data(x, y)
    return line,


anim = FuncAnimmap(animate, frames=M, interval=60, blit=True)
anim.save("heat_odesolver.gif")
