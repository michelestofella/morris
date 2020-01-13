# MORRIS-LECAR MODEL 

The goal of the repository is to perform the bifurcation analysis of the Morris Lecar model in order to discriminate neurons' excitability. 

The Morris Lecar model here used is described in the main paper (file `paper.pdf`) and the fixed parameters are defined in the file `MorrisLecar/fixed_parameters.py`. They are set as written in the main paper and can be modified to study different systems. Bifurcation analysis is performed by varying parameters `I_app` and `v_ca` that are not defined in `fixed_parameters.py` since they have to be parsed to the scripts. 

In order to properly perform bifurcation analysis, the Newton algorithm to find the zeros of a function is needed and so is an algorithm that integrates differential equations: in this case, the 4th order Runge Kutta algortihm was implemented. 
The repository can be thus divided into three parts:
1) Implementation of the Newton algorithm (folder `newton`)
2) Implementation of the Runge Kutta algorithm (folder `RungeKutta`)
3) Bifurcation analysis (folder `MorrisLecar`)

Both the Newton and the Runge Kutta algorithms have been tested via the libraries `pytest` and `hypothesis` and the tests are contained in the folder of the algorithm. Further description on the scripts and on the tests can be found in the folders `RungeKutta` for what concerns the Runge Kutta algorithm and `newton` for the Newton algorithm.

To reproduce figures like the ones in the main paper (here in folder `Images`), have a look at the folder `example` that shows how to obtain such figures under a particular configuration. 

All the scripts are written in **python**.
