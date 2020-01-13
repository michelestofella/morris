# RungeKutta

The present folder contains the implementation of the 4th order runge kutta algorithm to solve a system of n differential equations. The algorithm is contained in the file rk4_system.py.

The algorithm has been tested with several tests (see below) and the libraries pytest and hypothesis were used to develop testing strategies. 

A property of the algorithm, namely the proportionality of the global truncation error with respect to the 4th power of the integration time step, is implemented in the folder `error_analysis`.

A brief description of the scripts of the present folder is given below.

## rk4_system.py

The file **rk4_system.py** contains the 4th order Runge Kutta algorithm to solve a set of differential equations written in the form

``` 
dy1/dt = f1(t,y1,...,yn)
...
dyn/dt = fn(t,y1,...,yn)
```

To call the function, use the following command line:

`t, y = RK4_system(f, dt, y0, t0, Nstep)`

As input parameters, the function needs:
* f: the list of functions `[f1,...,fn]` in the left members of the system of differential equations. **Attention**: all the functions must be defined with variable t as first input variable;
* dt: the integration time step;
* y0: the list of initial conditions on the set of variables y1,...,yn;
* t0: initial condition on time;
* Nstep: number of iterations to be performed.

As output, the function returns:
* t: an Nstep-dimensional vector containing the time instants at which the solution has been calculated;
* y: a multidimensional array (size: `Nstep X Neq`, where `Neq` is the number of equations in the system) containing the solution of the system.

### rk4_system_test.py

The file **rk4_system_test.py** contains the tests of the 4th order Runge Kutta algorithm implemented in **rk4_system.py**. To perform the test, go to the folder RungeKutta and digit the following command line:

`pytest rk4_system_test.py`

We give here a description of the tests performed.

* `test_monodimensional_constant_case` integrates the simple differential equation `dy/dt=0` that has a constant value depending on the initial condition as solution. The test checks that the solution is constant. The strategy here applied is to see if this is true considering several initial conditions. Even the total number of steps is changed within the testing strategy. 
* `test_monodimensional_linear_case` integrates the simple differential equation `dy/dt=constant` that has as solution a line. The test checks if the solution is the expected line considering different initial guesses and different different total number of steps. 
* `test_bidimensional_independent_case` checks if the correct solution is reached while considering a bidimensional system formed by two differential equations, the first depending only on variable x, the second only on variable y. The strategy is to vary the initial conditions on both x and y and the total number of steps to be performed. 
* `test_bidimensional_dependent_case` integrates a simple system of two differential equations and checks that the correct solution is reached, varying the total number of steps and the initial conditions on both variables. 
* `test_reversibility` checks the reversibility of the algorithm in a linear monodimensional case. The algorithm is applied to go from time `0` to time `Nstep*dt`; then the time step is set as negative `-dt` and the algorithm is applied backward. The test checks if the initial value of the forward algorithm and the final value of the backward algorithm are equal. 

## rk4_error.py

This file is contained in the folder `error_analysis`.
To verify the correct implementation of the algorithm in **rk4.py**, we analyze the global truncation error for a simple model, namely the *passive membrane model* (see Ingalls (2013), *Mathamatical Modeling in Systems Biology. An Introduction*, MIT Press),

`dy/dt = g/C*(E-y)`

where g, C and E are model parameters. For such a simple model, an analitical solution exists:

`y(t) = E - np.exp(-g/C*t)*(E-v0)`

Thus we can compare the numerical and analitical solution. 

In the script, first we give an example of the integration. Then we evaluate the global truncation error. We set a total time for the integration and we apply the algorithm with different time steps.
We then evaluate the distance between the numerical solution at the last step and the analitical solution calculated at the last time step. 
We obtain that the global truncation error goes as the 4th power of the time step, as known from literature. 

![alt text](https://github.com/michelestofella/morris/blob/master/Images/rungekutta_error.png)

