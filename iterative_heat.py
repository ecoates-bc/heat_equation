# Solve the 1D heat equation with forward difference
# This is prone to instability!
# following this guide: http://people.uncw.edu/hermanr/pde1/NumHeatEqn.pdf

import numpy as np
import math
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

L = 1  # length of a rod
T = 0.1  # total time
k = 1  # heat constant
N = 20  # samples on the rod
M = 50  # timesteps

dx = L / N  # distance between samples on the rod
dt = T / M  # timestep

alpha = k * dt / dx**2

plane = np.zeros((M, N + 1))


def position_at(i):
    return i * dx


# initial condition
for i in range(N+1):
    plane[0, i] = math.sin(2 * math.pi * position_at(i))

# compute next timesteps
for m in range(1, M):
    for n in range(1, N):
        plane[m, n] = (
            plane[m - 1, n]
            + alpha
            * (plane[m - 1, n + 1] - 2 * plane[m - 1, n] + plane[m - 1, n - 1])
        )


# plot an animation
fig = plt.figure()
ax = plt.axes(xlim=(0,L), ylim=(-1,1))
plt.style.use("seaborn-pastel")
line, = ax.plot([], [], lw=3)


def animate(i):
    x = np.linspace(0, L, N+1)
    y = plane[i,:]
    line.set_data(x, y)
    return line,


anim = FuncAnimation(fig, animate, frames=M, interval=60, blit=True)
anim.save("heat_iterative.gif")
