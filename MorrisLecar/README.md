# MorrisLecar

This folder contains the automated scripts used to integrate the Morris Lecar model (`integrate.py`), to perform bifurcation analysis (`bifurcation_analysis.py`) and to draw the frequency plot (`frequency_plot.py`). For a description of the model, see [wikipedia](https://en.wikipedia.org/wiki/Morris–Lecar_model). For an accurate derivation of the model, see Ingalls (2013), *Mathematical Modeling in Systems Biology. An Introduction*, MIT Press. The parameters of the model are defined in the file `fixed_parameters.py`.

## fixed_parameters.py

In this file the needed libraries are imported, together with the 4th order Runge Kutta method (implemented in `morris\rungekutta\rk4_system.py`) and the bidimensional Newton algorithm (implemented in `morris\newton\newton2.py`). 

Furthermore, the model parameters are here set following, as in the main paper, Liu (2014), *Bifurcation Analysis of a Morris-Lecar neuruon model*, Biol Cybern 108;75-84. 

The parameters over which bifurcation analysis is to be implemented, namely `v_ca` and `I_app` are not defined in this file.

## integrate.py

This file contains the integration of the Morris Lecar model through the Runge Kutta algorithm (imported by `morris\rungekutta\rk4_system.py`) and a visualization of the signal voltage in time and of the phase space is shown. 
We give here a description of the script, following the block of code.

To run the integration, digit `python integrate.py` followed by the following parameters:
* --v_ca  the parameter of the morris lecar model that discriminates different classes of neurons.
* --I_app the external current applied to the model 
* --dt    integration time step
* --Nstep number of integration steps to be performed
* --v0    initial condition on the voltage
* --w0    initial condition of the fraction of opened channels
* --out   name of the generated figure

After parsing the parameters, we introduce the fixed parameters of the model as described by Liu (2014), importing them from `fixed_parameters.py`. 

Then, the Runge Kutta algorithm is applied following the parameters given as input.

Finally, a plot shows the evolution of the voltage through time and the trajectories in the phase space. The plot is saved as png file with the name given as input. 

## bifurcation_analysis.py

This file contains the bifurcation analysis of the model. In order to perform bifurcation analysis, we follow the following scheme:
* We set the value of `I_app`
* We apply the two dimensional Newton algorithm to find the zeros of the Morris Lecar model
* For each zero found, we calculate the eigenvectors of the Jacobian matrix to establish if the zero the point is stable or unstable.
* We change the value of `I_app`.

To run the script, digit the command line from the present directory `python bifurcation_analysis.py` followed by the parameters:
* --v_ca  the parameter of the morris lecar model that discriminates different classes of neurons.
* --Imin  minimum value of applied current analysed
* --Imax  maximum value of applied current analysed 
* --v0min minimum value of initial condition on voltage
* --v0max maximum value of initial condition on voltage
* --out   name of the generated figure

We describe here the script following the blocks of code. After parsing the parameters from the command line, the model parameters are imported from `fixed_parameters.py`, together with the bidimensional Newton algorithm (imported from `morris\newton\newton2.py`) and the Runge Kutta algorithm (imported from `morris\rungekutta\rk4_system.py`).

Then, we set the range of values of `I_app` and the set of initial conditions `v0` to perform bifurcation analysis; depending on the parsed parameters, we have that `I_app` lies in the interval `[Imin,Imax]` while the initial condition `v0` in `[v0min,v0max]`. 

Now we can start bifurcation analysis: for each value of `I_app`, we find all the zeros of the model by exploiting several initial guess values and applying the bidimensional Newton algorithm. 
All the zeros found are appended to two lists `v_zeros` and `w_zeros`.

In the next block, stability is discriminated. We calculate the eigenvalues of the Jacobian matrix for each zero previously found. If all the eigenvalues are negative, than the point is stable, otherwise it is unstable. 
Two lists, namely `stable` and `unstable` are generated and contain respectively the stable and unstable points. Each element of the list is itself a list containing `I_app, v_zeros, w_zeros, stability`, where the last element is a string containing *stable* or *unstable*.

Before plotting the data, we reorganize them into pandas dataframes so that they are more tidy.

Finally, we plot the bifurcation diagram: for each value of `I_app`, we plot the stable (solid line) and the unstable (dashed line) points. 
The figure is saved as png file with a name given by the input parameter `--out`.

**Attention.** The algorithm often shows errors while running. This happens because some initial guesses or some values of applied current can lead the Newton algorithm not to converge because the determinant of the jacobian in that particular point is zero or because the function has no zeros with that particular value of applied current. See the implementation of the Newton algorithm for furhther information (folder `morris/newton`). 
Even if such exceptions are shown to let the user be aware of the problem, the points associated to these exceptions are not considered in the final plot. Summarizing, if an exception arises from the Newton algorithm, the point is skipped. 

## frequency_plot.py

This file contains the script to reproduce the frequency plot. 
To calculate the frequency of the generated output signal, we use the Scipy function `find_peaks`.

To run the script, digit the command `python frequency_plot.py` followed by the parameters below:
* --v_ca  the parameter of the morris lecar model that discriminates different classes of neurons.
* --dt    integration time step
* --Nstep number of integration steps to be performed
* --v0    initial condition on the voltage
* --w0    initial condition of the fraction of opened channels
* --Imin  minimum value of applied current analysed
* --Imax  maximum value of applied current analysed 
* --out   name of the generated figure

After parsing the parameters, we need to import the parameters and the Runge Kutta algorithm from `fixed_parameters.py` and we have to define the Morris Lecar model for the Runge Kutta integration.

We then have to set the values of `I_app` through which we want to calculate the signal. The parameters for the integration are defined by the parsed parameters and the Runge Kutta algorithm is applied.
At the end of each application of the algorithm, the generated signal is analyzed. If the action potential has been generated, then the generated signal is periodic. Using the function `find_peaks`, we detect the maxima of the signal, we take the mean value of the periods of the signal and we take as frequency the inverse of such a period. 

**Attention.** The action potential is generated above a threshold. For class1 neurons (see `paper.pdf` for further details), the output frequency should be arbitrarily small. Above the threshold, the higher the `I_app`, the higher the frequency. To reach smaller frequencies, the system should be integrated for an appropriate time interval, otherwise the second peak cannot be reached. 

Finally, the frequency plot is shown: for each `I_app` value, the calculated frequency is plotted. The figure is saved as a png file with the name given with the input parameter `--out`. 
