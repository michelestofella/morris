# =============================================================================
# ERROR ANALYSIS OF THE RUNGE-KUTTA ALGORITHM
# 
# The program performs the analysis of the global truncation error
# of the Runge-Kutta algorithm for the specific (and simple) case of the 
# 'Passive Cell Membrane Model'. 
#
# Algorithm source:     https://rosettacode.org/wiki/Runge-Kutta_method
# Algorithm tests:      !pytest rk4.py
# 
# Source for the model: INGALLS, 'Mathematical Modeling in Systems Biology',
#                       The MIT Press (2013)
# 
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt    

""" Runge-Kutta algorithm """
def RK4(f, dt, y0, t0, Nstep):
    y = np.zeros(Nstep+1); t = np.zeros(Nstep+1)
    t[0] = t0; y[0] = y0
    for i in range(0,Nstep):
        dy1 = dt*f(t, y)
        dy2 = dt*f(t+0.5*dt, y+0.5*dy1)
        dy3 = dt*f(t+0.5*dt, y+0.5*dy2)
        dy4 = dt*f(t+dt,     y+dy3)
        
        y[i+1] = y[i] + (dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])/6
        t[i+1] = t[i] + dt
    return t,y 

# %%

""" PASSIVE CELL MEMBRANE MODEL """
g = 0.0144; C = 0.98; E = -93.6; 
dt = 0.1; v0 = -60.0; t0 = 0.0; Nstep = 5000

def single_ion(t,y):
    return g/C*(E-y)  

def solution(t):
    return E - np.exp(-g/C*t)*(E-v0)

time, voltage = RK4(single_ion, dt, y0=v0, t0=t0, Nstep=Nstep) 
    
plt.plot(time, voltage, label='Numerical', c='black')
plt.plot(time, solution(time), 'r--', label='Analitical')  
plt.xlabel('Time [ms]', fontsize=15)
plt.ylabel(' Membrane Voltage [mV]', fontsize=15) 
plt.legend()     

# %%

""" ERROR ANALYSIS OF THE RUNGE KUTTA ALGORITHM
    The global truncation error should go ad dt^4 """
""" We calculate the Global Truncation Error """
total_time = 100.0
time_steps = np.linspace(0.01,1.0,20)
Nsteps = np.round(total_time/time_steps)
Nsteps = Nsteps.astype(np.int64)
error = np.zeros(len(time_steps))

for i in range(0, len(time_steps)):
    time, voltage = RK4(single_ion, dt=time_steps[i],
                        y0=v0, t0=t0, Nstep=Nsteps[i])
    error[i] = np.abs(solution(time[-1])-voltage[-1])  
# %%
""" Linear Fit of Error vs dt^4 """
dt4 = np.power(time_steps, 4)    
coef = np.polyfit(dt4, error, 1)
poly = np.poly1d(coef)
# %%
""" Visualization of the results """
xx = np.linspace(min(dt4),max(dt4))
plt.scatter(dt4, error, c='black')
plt.plot(xx, poly(xx), 'r--')
plt.ylim(min(error),max(error))
plt.xlabel('$\Delta$t$^4$', fontsize=15)
plt.ylabel('Error', fontsize=15)

# %%