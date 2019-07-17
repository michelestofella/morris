# MorrisLecar

This folder contains the scripts used to integrate the Morris Lecar model. For a description of the model, see [wikipedia](https://en.wikipedia.org/wiki/Morris–Lecar_model). For an accurate derivation of the model, see Ingalls (2013), *Mathematical Modeling in Systems Biology. An Introduction*, MIT Press. 

## morris_setup.py

In the file **morris_setup.py** the needed libraries are imported, together with the 4th order Runge Kutta method (implemented in `morris\rungekutta\rk4_system.py`) and the bidimensional Newton algorithm (implemented in `morris\newton\newton2.py`). 

Furthermore, the model parameters are here set following Liu (2014), *Bifurcation Analysis of a Morris-Lecar neuruon model*, Biol Cybern 108;75-84. 

The parameter changed to perform bifurcation analysis in the paper is the parameter named `v_ca`.

## morris.py

The file **morris.py** contains the integration of the Morris Lecar model through the Runge Kutta algorithm (imported by `morris\rungekutta\rk4_system.py`) and a visualization of the signal voltage in time and of the phase space is shown. 

First of all, we introduce the model as described by Liu (2014), importing the parameters from the **morris_setup.py** file. 

Then, the Runge Kutta algorithm is applied for two different initial conditions.

Finally, a plot shows the evolution of the voltage through time and the trajectories in the phase space. An examle is given by the following image.

![alt text](https://github.com/michelestofella/morris/blob/master/Images/class2(i%3D80%2Cv0%3D-25%2C-50%2Cw0%3D0%2C01).png)

## morris2.py

The file **morris2.py** contains the bifurcation analysis of the model. In order to perform bifurcation analysis, we follow the following scheme:
* We set the value of `I_app`
* We apply the two dimensional Newton algorithm to find the zeros of the Morris Lecar model
* For each zero found, we calculate the eigenvectors of the Jacobian matrix to establish if the zero the point is stable or unstable.
* We change the value of `I_app`, covering the range 0-100

First of all, the model parameters are imported from **morris_setup.py**, together with the bidimensional Newton algorithm (imported from `morris\newton\newton2.py`) and the Runge Kutta algorithm (imported from `morris\rungekutta\rk4_system.py`).
The Morris Lecar model is then described in two versions:
* the model for the Runge Kutta algorithm needs an explicit dependence on time in the input variables
* the model for the Newton algorithm must avoid the time dependence and requires the Jacobian matrix

Then, we set the range of values of `I_app` and the set of initial conditions `v0` to perform bifurcation analysis. For example,
```
I_app_values = np.linspace(0,100,201)    
v0_values = np.linspace(-80,40,61)
```
Now we can start bifurcation analysis: for each value of `I_app`, we find all the zeros of the model by exploiting several initial guess values and applying the bidimensional Newton algorithm. 
All the zeros found are appended to two lists `v_zeros` and `w_zeros`.

In the next block, stability is discriminated. We calculate the eigenvalues of the Jacobian matrix for each zero found previously. If all the eigenvalues are negative, than the point is stable, otherwise it is unstable. 
Two lists, namely `stable` and `unstable` are generated and contain respectively the stable and unstable points. Each element of the list is itself a list containing `I_app, v_zeros, w_zeros, stability`, where the last element is a string containing *stable* or *unstable*.

Since the model shows the birth of a limit cycle above a certain threshold, in particular when a stable point becomes unstable through a SNIC bifurcation or a Hopf bifurcation, we integrate the model using the Runge Kutta algorithm in order to show the maximum and minimum values taken by the voltage variable `v` in the limit cycle. 
**Attention**: to get reasonable results, the system should be integrated for an appropriate time interval, so that the trajectory in the phase space is able to cover all the values of the limit cycle. 

Before plotting the data, we reorganize them into pandas dataframes so that they are more ordered.

Finally, we plot the bifurcation diagram: for each value of `I_app`, we plot the stable (solid line) and the unstable (dashed line) points. The eventual bifurcation is indicated by a vertical dotted line, while the maximum and minimum values assumed by the voltage variable `v` in the limit cycle are solid lines labeled with `Vmax` and `Vmin`.
An example is here shown.

![alt text](https://github.com/michelestofella/morris/blob/master/Images/class2_bif.png)

## morris3.py

The file **morris3.py** contains the script to reproduce the frequency plot. 

To calculate the frequency of the generated output signal, we use the Scipy function `find_peaks`.

First of all, we need to import the parameters and the Runge Kutta algorithm from **morris_setup.py** and we have to define the Morris Lecar model for the Runge Kutta integration.

We then have to set the values of `I_app` through which we want to calculate the signal. The parameters for the integration are defined and the Runge Kutta algorithm is applied.
At the end of each application of the algorithm, the generated signal is analyzed. If the action potential has been generated, then the generated signal is periodic. Using the function `find_peaks`, we detect the maxima of the signal, we take the mean value of the periods of the signal and we take as frequency the inverse of such a period. 

**Attention**: the action potential is generated above the threshold. For class1 neurons (see the paper for further details), the output frequency should be very small. Above the threshold, the higher the `I_app`, the higher the frequency. To reach smaller frequencies, the system should be integrated for an appropriate time interval, otherwise the second peak is not reached. 

Finally, the frequency plot is shown: for each `I_app` value, the calculated frequency is plotted. An example is shown in the following image. 

![alt text](https://github.com/michelestofella/morris/blob/master/Images/class1_freq.png)

