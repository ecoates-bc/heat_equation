# Solving the 2D heat equation with dune-fem

import matplotlib
matplotlib.rc("image", cmap="jet")

import numpy as np
from scipy.sparse.linalg import spsolve as solver

from dune.grid import structuredGrid, GridFunction
from dune.fem import assemble
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin as solutionScheme
from dune.ufl import DirichletBC, Constant
from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    dx,
    grad,
    dot,
    exp,
)

T = 10.0  # final time
num_steps = 50  # number of steps
dt = T / num_steps  # timestep size

# define the 2D plate
gridView = structuredGrid([-2, -2], [2, 2], [30, 30])
space = lagrange(gridView, order=2)
# boundary conditions
bc = DirichletBC(space, Constant(0))

# Initial condition: gaussian hill
x = SpatialCoordinate(space)
a = 5
initial = exp(-1*a*(x[0]**2) - a*(x[1]**2))

u_h = space.interpolate(initial, name="u_h")
u_h_n = u_h.copy(name="previous")


# define trial and test functions
u = TrialFunction(space)
v = TestFunction(space)
f = Constant(0)

t = Constant(0, name="t")

# define weak form
# TODO: is there a dune equivalent to fenics' a, L = lhs(F), rhs(F) ?
a = u*v*dx + dt*dot(grad(u), grad(v))*dx
L = (u_h + dt * f) * v * dx

scheme = solutionScheme([a==L, bc], solver="cg")

for n in range(num_steps):
    t += dt
    u_h_n.assign(u_h)
    y = u_h.as_numpy
    y[:] = solver(A, b)
    u_h.assign(y)
    u_h.plot()