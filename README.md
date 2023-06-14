# Some experiments with PDE solving

Here's some initial work I did with FDM/FEM and the Heat equation. It's a mix of tutorials I followed and some attempts at writing my own code.

## Finite Difference
I started with following a tutorial to make `iterative_heat.py` and write a simple iterative FD solution to the heat equation. The result is in `heat_iterative.gif` - it suffers from numerical stability when the function approaches zero.

The next thing I tried was turning the equation into a system of ODEs in `heat_system.py`. The results are in `heat_odesolver.gif`. Looks like the system is a lot more, if not completely, stable.

## FEniCSx
The first FEA library I used was FEniCSx. I was able to get a simple 2D solution to the heat equation working, in `heat_equation_fenics.ipynb`.

## dune-fem
The next library I tried was DUNE - I followed some tutorials, and got close to getting a 2D heat equation solution working, but I'm not confident that my approach worked (mainly because the visualizations at the end aren't actually moving).

## FEniCSx vs DUNE
I found that FEniCSx was overall easier to program with and faster than DUNE. However, it was trickier to get working on my computer. I ended up having to pull a docker image and work in a jupyter environment - I have a feeling that FEniCSx would be quite tricky to get up and running on Sockeye, for example. DUNE was much more straightforward to get up and running, but is slower (from my short time using it so far.) It looks like both have parallel support, but I have a feeling that DUNE will be more straightforward to run in parallel.