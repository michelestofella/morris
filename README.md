# MORRIS 

Analysis of Bifurcation in the Morris Lecar Model.

Project for 'Models and Numerical Methods for Applied Physics'.

## RungeKutta

The folder contains the implementation of the Runge-Kutta algorithm and the error analysis of the algorithm applied to a simple model, 
namely the passive cell membrane potential.

### rk4.py

The program RK4.PY contains the implementation of the Runge-Kutta algorithm (source:https://rosettacode.org/wiki/Runge-Kutta_method).
The algorithm aims to solve differential equations written as

  dy/dt = f(t,y)
  
To call the algorithm:
  
  RK4 (f, dt, y0, t0, Nstep)

Parameters
- f(t,y)  is a function that may depend both on y and t;
- dt      is the integration time step;
- y0      is the initial condition;
- t0      is the initial time;
- Nstep   is the number of steps for which the algorithm must be repeated

  
The program also contains two tests:

  1) the first test applies the algorithm forward and backward in time and sees if the final value is equal to the initial one 
  within a certain tolerance on precision governed by the parameter tol (set at tol=1e-5);
  
  2) the second test applies the algorith three times: the first one to go from A to B, then to go from B to C and then to go from A to C. 
  The test shows that the final point of the algorithm applied A -> B and B -> C is equal to the algorithm applied directly A -> C 
  within a certain tolerance on precision governed by the parameter tol (set at tol=1e-5).  
