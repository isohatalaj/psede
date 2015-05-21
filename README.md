# psede

A library for solving differential equations using pseudospectral
methods.

The design goal has been to create a set of routines that allow for a
quick and easy solution of ordinary (and later, partial) differential
equations for which the usual local iterative methods (Runge-Kutta,
predictor-correctors, and friends) are not well suited. Such problems
are encountered when eg. solving for unstable solutions, solutions
that live on separatrices, singular equations, etc. 

As general purpose routines, they are probably not the best possible
choice for any given problem. The idea here is to rather provide
something "good enough" that works, can be set up quickly, and can
provide a basis for more refined approaches.

This library was motivated by continuous time optimal control problems
in economics, specifically solving the Hamilton-Jacobi-Bellman
equation, a nonlinear ODE/PDE with a variety of boundary conditions.




