# PDE Solution to Krusell-Smith (1998)

This is a Matlab program which uses **Finite-Difference Method** to solve the Krusell Smith (1998) model. Krusell-Smith (1998) proposed a 
multiple-agent model with uninsurable idiosyncratic income shocks and aggregate shock. The algorithm first use a continuous-time version 
of Krusell-Smith model, and solve it using the finite-difference method.

We first write down the model:

```math
V(a_t, z_t) = \max_{c_t} \int_0^\infty e^{-\rho t} u(c_t) dt
```
subject to budget constraint:

```math
\dot{a}_t = r_t a_t + y_t - c_t
```

**Stationary Equilibrium** is defined to be the equilibrium in which the distribution of the population is stationary. Hence the interest rate
is constant, and we can write the above problem into a Hamiltonian-Jaccobi Bellman Equation:

```math 
\rho V(a) = \max_{c} \{ u(c) + \partial_a V(a, z)(ra + y - c) + \lambda(V(a,z') - V(a,z))\}
```
