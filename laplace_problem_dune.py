# Solving a laplace problem with dune-fem
# taken from https://www.dune-project.org/sphinx/content/sphinx/dune-fem/dune-fempy_nb.html

import matplotlib
matplotlib.rc("image", cmap="jet")

import numpy as np
from scipy.sparse.linalg import spsolve as solver

from dune.grid import structuredGrid
from dune.fem import assemble
from dune.fem.space import lagrange
from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    dx,
    grad,
    inner,
    dot,
    sin,
    cos,
    pi
)


gridView = structuredGrid([0,0], [1, 1], [20,20])
space = lagrange(gridView, order=1)
u_h = space.interpolate(0, name="u_h")

x = SpatialCoordinate(space)
u = TrialFunction(space)
v = TestFunction(space)

f = (8*pi**2 + 1) * cos(2*pi*x[0])*cos(2*pi*x[1])
a = (inner(grad(u), grad(v)) + u*v) * dx
l = f*v * dx

mat, rhs = assemble(a==l)

A = mat.as_numpy
b = rhs.as_numpy
y = u_h.as_numpy
y[:] = solver(A, b)

u_h.plot()