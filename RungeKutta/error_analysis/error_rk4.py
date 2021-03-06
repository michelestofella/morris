# =============================================================================
#
# ERROR ANALYSIS OF THE RUNGE-KUTTA ALGORITHM
# 
# The code performs the analysis of the global truncation error
# of the Runge-Kutta algorithm for the specific (and simple) case of the 
# 'Passive Cell Membrane Model'. 
# 
# Source for the model: INGALLS, 'Mathematical Modeling in Systems Biology',
#                       The MIT Press (2013)
# 
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt    
from rk4_system import RK4_system

# %% 

""" Passive cell membrane model """
def single_ion(t,y):
    return g/C*(E-y)  
""" Analitical solution """
def solution(t):
    return E - np.exp(-g/C*t)*(E-v0)

""" Model parameters """
g = 0.0144; C = 0.98; E = -93.6; 
""" Integration parameters """
dt = 0.1; v0 = -60.0; t0 = 0.0; Nstep = 5000

time, voltage = RK4_system([single_ion], dt, [v0], t0=t0, Nstep=Nstep) 
    
plt.plot(time, voltage[0], label='Numerical', c='black')
plt.plot(time, solution(time), 'r--', label='Analitical')  
plt.xlabel('Time [ms]', fontsize=15)
plt.ylabel(' Membrane Voltage [mV]', fontsize=15) 
plt.legend()     

# %%

""" ERROR ANALYSIS OF THE RUNGE KUTTA ALGORITHM
The global truncation error should go ad dt^4 """

total_time = 100.0
time_steps = np.linspace(0.01,1.0,20)
Nsteps = np.round(total_time/time_steps)
Nsteps = Nsteps.astype(np.int64)
error = np.zeros(len(time_steps))

""" We calculate the global truncation error for each time step """
for i in range(0, len(time_steps)):
    time, voltage = RK4_system([single_ion], dt=time_steps[i],
                        y0=[v0], t0=t0, Nstep=Nsteps[i])
    error[i] = np.abs(solution(time[-1])-voltage[0][-1])  

""" Linear Fit of Error vs dt^4 """
dt4 = np.power(time_steps, 4)    
coef = np.polyfit(dt4, error, 1)
poly = np.poly1d(coef)

""" Visualization of the results """
xx = np.linspace(min(dt4),max(dt4))
plt.scatter(dt4, error, c='black')
plt.plot(xx, poly(xx), 'r--')
plt.ylim(-0.5e-9,5e-9)
plt.xlabel('$\Delta$t$^4$', fontsize=25)
plt.ylabel('Error', fontsize=25)
plt.tick_params(axis='both', labelsize=15)
plt.grid(True)

# %%