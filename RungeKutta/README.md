# RungeKutta

## rk4.py

The file **rk4.py** contains the 4th order Runge Kutta algorithm to solve one dimensional differential equation written in the form

`dy/dt = f(t,y)`

For a description of the algorithm see  [wikipedia](https://en.wikipedia.org/wiki/Runge–Kutta_methods) or [rosettacode](https://rosettacode.org/wiki/Runge-Kutta_method).

To call the function, use the following command:

`t, y = RK4(f, dt, y0, t0, Nstep)`

As input parameters, the function needs:
* f: the function f(t,y) in the left side of the differential equation. **Attention**: the function must be defined with variable t as first input parameter;
* dt: the integration step;
* y0: initial condition of the y variable;
* t0: initial condition on time variable;
* Nstep: number of iterations that must be performed;

As output the function returns:
* t: an Nstep-dimensional vector containing the time istants at which the solution has been calculated; the difference between two adjacent elements of this vector is dt; 
* y: an Nstep-dimensional vector containing the solution of the differential equation.

### rk4_test.py

The file **rk4_test.py** contains the tests of the 4th order Runge Kutta algorithm implemented in **rk4.py**. To perform the test, you should go to the folder RungeKutta and digit the following command line in the iPython console:

`!pytest rk4_test.py`

We here give a description of each test. The file contains only one test. Further testing of the algorithm is available in the **rk4_error.py** file where error analysis is performed. 

* `test_one` applies the algorithm to the differential equation `dy/dt=0` and checks if the last element of the returned solution is equal to the initial value: `assert y[-1]==y0`

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

The file **rk4_system_test.py** contains the tests of the 4th order Runge Kutta algorithm implemented in **rk4_system.py**. To perform the test, go to the folder RungeKutta and digit the following command line in the iPython console:

`!pytest rk4_system_test.py`

We give here a description of the tests performed.

* `test_1` applies the algorithm to the differential equation `dy/dt=0` and checks if the last element of the solution is equal to the initial value: `assert y[0][-1]==y0[0]`.
* `test_2` applies the algorithm to the differential equation `dy/dt=1` with initial condition `y0=0` and checks if the last element of the t array and the solution array is equal: `assert y[0][-1]==t[-1]`.
* `test_3` applies the algorithm to the system of two differential equation `dy1/dt=0` and `dy2/dt=1` and checks if the last element of the solution of the first differential equation is equal to the initial value `y[0][-1]=y0[0]` and if the final value of the solution of the second differential equation is equal to the final value of the time array `y[1][-1]==t[-1]`.
* `test_4` applies the algorithm to the system of two differential equations `dy1/dt=-y1` and `dy2/dt=-y2` and checks if the last element of the two solutions is equal: `assert y[0][-1]==y[1][-1]`.
* `test_5` tests the reversibility of the algorithm in the linear case `dy/dt=-2`. The algortihm is applied forward for `Nstep=1000`; then the last value is set as initial condition for the backward algortihm (i.e. the algorithm applied with inverse time step `dt=-dt`). The test checks if the last value is equal to the initial condition. 


## rk4_error.py

To verify the correct implementation of the algorithm in **rk4.py**, we analyze the global truncation error for a simple model, namely the *passive membrane model* (see Ingalls (2013), *Mathamatical Modeling in Systems Biology. An Introduction*, MIT Press),

`dy/dt = g/C*(E-y)`

where g, C and E are model parameters. For such a simple model, an analitical solution exists:

`y(t) = E - np.exp(-g/C*t)*(E-v0)`

Thus we can compare the numerical and analitical solution. 

In the script, first we give an example of the integration. Then we evaluate the global truncation error. We set a total time for the integration and we apply the algorithm with different time steps.
We then evaluate the distance between the numerical solution at the last step and the analitical solution calculated at the last time step. 
We obtain that the global truncation error goes as the 4th power of the time step, as known from literature. 

![alt text](https://github.com/michelestofella/morris/blob/master/Images/rungekutta_error.png)

