# Bifurcation Analysis of Morris Lecar Model

The goal of the repository is to perform the bifurcation analysis of the Morris Lecar model in order to discriminate neurons' excitability. The results of the analysis are discussed in the file `paper.pdf`.

All the scripts are written in **python** (version 3.7.3). Libraries needed: `argparse`, `hypothesis`, `matplotlib`, `numpy`, `pandas`, `pytest`, `scipy` and `sys`.

## Bifurcation Analysis

In dynamical systems, a bifurcation occurs when a small smooth change made to a parameter (the bifurcation parameter) of a system causes a sudden “qualitative” or topological change in its behaviour. For further information about bifurcation theory, look at [wikipedia](https://en.wikipedia.org/wiki/Bifurcation_theory). 

In order to properly perform bifurcation analysis, we need an algorithm to find the zeros of a function and an algorithm that integrates differential equations: in this case the [Newton algorithm](https://en.wikipedia.org/wiki/Newton%27s_method) the 4th order [Runge Kutta algorithm](https://en.wikipedia.org/wiki/Runge–Kutta_methods) are implemented from scratch. 
The repository can be thus divided into three parts:
1) Implementation of the Newton algorithm (folder `newton`); further documentation [here](https://github.com/michelestofella/morris/blob/master/newton/README.md).
2) Implementation of the Runge Kutta algorithm (folder `RungeKutta`); further documentation [here](https://github.com/michelestofella/morris/blob/master/RungeKutta/README.md).
3) Automated bifurcation analysis (folder `MorrisLecar`); further documentation [here](https://github.com/michelestofella/morris/blob/master/MorrisLecar/README.md).

Both the Newton and the Runge Kutta algorithms have been tested via the libraries `pytest` and `hypothesis`. Further description on the scripts and on the tests can be found in the folders `RungeKutta` for what concerns the Runge Kutta algorithm and `newton` for the Newton algorithm.

## Morris Lecar Model

The [Morris Lecar model](https://en.wikipedia.org/wiki/Morris–Lecar_model) here used is described in the main paper (file `paper.pdf`). It describes the time evolution of the voltage in a membrane domined by two ionic currents: an inward calcium current and an outward potassium current. Experimental studies show that calcium channels relax to steady state much more quickly than the potassium channels, thus the model can be described by two only variables, namely the voltage `V` and the fraction of opened potassium channels `w`.

Many parameters have to be defined. The parameters over which bifurcation analysis is performed:

* `I_app` is a current applied which describes any current injected into the cell by an experimenter or generated in response to signals from other cells;
* `v_ca` is a numerical parameter that determines the steady state of calcium channels;

and some fixed parameters:

* `C` is the conductance of the membrane; default: `C=2`
* `g_ca` and `g_k` are the maximal membrane conductance for calcium and potassium channels; default: `g_ca=20` and `g_k=20`
* `E_ca` and `E_k` are the Nernst potential of calcium and potassium channels; default: `E_ca=50` and `E_k=-100'
* `g_leak` and `E_leak` are the maximal membrane conductance and the Nernst potential of a leak current introduced to take into account the background activity of other (secundary) ion fluxes; default: `g_leak=2` and `E_leak=-70`
* `v_k`,`theta_ca` and `theta_k` are other numerical parameters that together with `v_ca` determine the steady state of the calcium and potassium channels; default `v_k=-10`, `theta_ca=18` and `theta_k=13`
* `phi_w` is a constant term introduced by [Liu](https://www.researchgate.net/publication/259770887_Bifurcation_analysis_of_a_Morris-Lecar_neuron_model); default: `phi_w=0.15`.

All the fixed parameters are defined in the file `MorrisLecar/fixed_parameters.py`. They are set as written in the main paper and can be modified to study different systems. Bifurcation analysis is performed by varying two parameters, namely `I_app` and `v_ca`, that are not defined in `fixed_parameters.py` since they have to be parsed to the scripts. 

## Example

The folder `example` shows how to properly run the automated scripts to integrate the model, to perform bifurcation analysis and to show the frequency plot. If the model is integrated via the Runge Kutta algorithm, the time evolution of generated signal is shown together with the phase space; if bifurcation analysis is performed, the bifurcation diagram is plotted showing stable and unstable points of the model; the frequency plot shows how the frequency of the generated signal depends on the current applied to the system `I_app` if the generated signal is periodic.
