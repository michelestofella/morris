# MORRIS-LECAR MODEL 

Analysis of Bifurcation in the Morris Lecar Model.

Project for 'Models and Numerical Methods for Applied Physics'.

## RungeKutta

The folder contains all the codes about the implementation of the 4th order Runge-Kutta algorithm used to integrate the system of differential equations that represents the Morris-Lecar model. 

In particular, you should find the following files
* **rk4.py**: it contains the Runge-Kutta algorithm to solve one differential equation;
* **rk4_error.py**: it contains the analysis of the global truncation error introduced by the algorithm;
* **rk4_system.py**: it contains the Runge-Kutta algorithm to solve a set of two differential equations.

### rk4.py

The program **rk4.py** contains the implementation of the 4th order Runge-Kutta algorithm to solve one differential equation (source: <https://rosettacode.org/wiki/Runge-Kutta_method>).
The algorithm aims to solve differential equations written as
> dy/dt = f(t,y)
  
To call the algorithm:

`t, y = RK4 (f, dt, y0, t0, Nstep)`

**Parameters:**
* **f(t,y)**: is a function that may depend both on y and t;
* **dt**: is the integration time step;
* **y0**: is the initial condition;
* **t0**: is the initial time;
* **Nstep**: is the number of steps for which the algorithm must be repeated

**Return Value:**
* **t**: it is a float64 array of size *Nstep*
* **y**: it is a float64 array of size *Nstep* that represents the solution of the differential equation
  
The program **rk4.py** also contains two tests:

The **first test** applies the algorithm forward and backward in time (in both cases the algorithm is repeated *Nstep* times) to a test function 
  
  ```
  import numpy as np
  f(t,y) = t*np.sqrt(y)
  ```
  and checks if the final value is equal to the initial one within a certain tolerance on precision governed by the parameter tol (set at tol=1e-6);
  
The **second test** applies the algorith three times: the first time the algorithm is applied to go from an initial time value A to B and then to go from B to C; then the algorithm is applied to go directly from A to C. The test function is again
  
  ```
  f(t,y) = t*np.sqrt(y)
  ```
  
The test checks if the final point of the algorithm applied A -> B and B -> C is equal to the algorithm applied directly A -> C within a certain tolerance on precision governed by the parameter tol (set at tol=1e-6).  

**Note on the tests**: the two tests correctly work because an appropriate integration step dt has been set. See the file **rk4_error.py** for further details on the dependance of the error of the algorithm with respect to the integration step.

### rk4_error.py

### rk4_system.py
