# MORRIS-LECAR MODEL 

Analysis of bifurcations in the Morris Lecar model for the classification of neurons' excitability.

Project for 'Models and Numerical Methods for Applied Physics'.

## Structure of the repository

* **Bifurcation** contains some simple examples of bifurcation analysis.
* **Images** contains the images of the article. 
* **MorrisLecar** contains the scripts for the final analysis of the Morris Lecar model.
* **RungeKutta** contains the implementation and the tests of the 4th order Runge Kutta algorithm.
* **newton** contains the implementation and the tests of the Newton algorithm.
* **README.md** contains a description of the whole repository.

### Bifurcation

* **bifurc.py** contains some simple examples of bifurcation analysis.

### MorrisLecar

* **morris.py** contains the numerical integration of the Morris Lecar model.
* **morris2.py** contains the script to perform bifurcation analysis, to built the frequency plot and to integrate the model to draw the images of the paper.

### RungeKutta

* **rk4.py** contains the Runge-Kutta algorithm to solve one differential equation and its tests.
* **error_rk4.py** contains the analysis of the global truncation error introduced by the algorithm.
* **rk4_system.py** contains the Runge-Kutta algorithm to solve a set of two differential equations and its tests.

### newton 
* **newton.py** contains the one-dimensionale Newton algorithm and its tests.
* **newton2.py** contains the two-dimensional Newton algorithm and its tests. 
