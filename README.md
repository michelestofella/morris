# MORRIS-LECAR MODEL 

Analysis of bifurcations in the Morris Lecar model for the classification of neurons' excitability.

Project for 'Models and Numerical Methods in Physics'.

Language: **python**.

IDE: **spyder**

## Structure of the repository

* **Images** contains the images of the article. 
* **MorrisLecar** contains the scripts for the final analysis of the Morris Lecar model.
* **RungeKutta** contains the implementation and the tests of the 4th order Runge Kutta algorithm.
* **newton** contains the implementation and the tests of the Newton algorithm.
* **README.md** contains a description of the whole repository.
* **paper.pdf** contains the final paper in pdf format.

**Note**: every folder where some scripts have been implemented contains a *README* file with detailed documentation.

### MorrisLecar

* **morris.py** contains the numerical integration of the Morris Lecar model.
* **morris2.py** contains the script to perform bifurcation analysis.
* **morris3.py** contains the script to built the frequency plot
* **morris_setup.py** contains the parameters of the model

### RungeKutta

* **rk4.py** contains the Runge-Kutta algorithm to solve one differential equation.
* **rk4_test.py** contains the test of the Runge-Kutta algorithm implemented in rk4.py.
* **error_rk4.py** contains the analysis of the global truncation error introduced by the algorithm.
* **rk4_system.py** contains the Runge-Kutta algorithm to solve a set of two differential equations.
* **rk4_system_test.py** contains the test of the Runge-Kutta algorithm implemented in rk4_system.py.


### newton 
* **newton.py** contains the one-dimensionale Newton algorithm.
* **newton_test.py** contains the tests of the Newton algorithm implemented in newton.py
* **newton2.py** contains the two-dimensional Newton algorithm. 
* **newton2_test.py** contains the tests of the Newton algorithm implemented in newton2.py

