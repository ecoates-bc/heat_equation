# Solvling the Forchheimer problem with Dune
# taken from https://www.dune-project.org/sphinx/content/sphinx/dune-fem/concepts_nb.html#A-Nonlinear-Time-Dependent-Problem

import matplotlib
import numpy as np

from dune.grid import structuredGrid
from dune.fem.space import lagrange
from dune.fem.scheme import galerkin
from dune.ufl import Constant
from ufl import (
    TestFunction,
    TrialFunction,
    SpatialCoordinate,
    FacetNormal,
    ds,
    div,
    dot,
    sqrt,
    exp,
    dx,
    inner,
    grad,
)

gridView = structuredGrid([0,0], [1,1], [20,20])
space = lagrange(gridView, order=2)
u_h = space.interpolate(0, name="u_h")
u_h_n = u_h.copy(name="previous")

x = SpatialCoordinate(space)
initial = 1/2*(x[0]**2+x[1]**2) - 1/2*(x[0]**3 - x[1]**3) + 1

dt = Constant(0, name="dt")
t = Constant(0, name="t")
u = TrialFunction(space)
v = TestFunction(space)

abs_du = lambda u: sqrt(inner(grad(u), grad(u)))
K = lambda u: 2/(1 + sqrt(1 + 4*abs_du(u)))
a = (dot((u - u_h_n)/dt, v) \
    + 0.5*dot(K(u)*grad(u), grad(v)) \
    + 0.5*dot(K(u_h_n)*grad(u_h_n), grad(v))) * dx

exact = lambda t: exp(-2*t)*(initial - 1) + 1
f = lambda s: -2*exp(-2*s)*(initial-1) - div(K(exact(s))*grad(exact(s)))
g = lambda s: K(exact(s))*grad(exact(s))
n = FacetNormal(space)
b = 0.5*(f(t)+f(t+dt))*v*dx + 0.5*dot(g(t)+g(t+dt),n)*v*ds

scheme = galerkin(a == b, solver="cg")

endTime = 0.25
exact_end = exact(endTime)

scheme.model.dt = 0.01
time = 0
u_h.interpolate(initial)
step = 0
while time < (endTime - 1e-6):
    scheme.model.t = time
    u_h_n.assign(u_h)
    scheme.solve(target=u_h)
    time += scheme.model.dt
    step += 1

    if step % 5 == 0:
        u_h.plot()
