# =============================================================================
# SECOND VERIONS OF RUNGE-KUTTA ALGORITHM
# This second version uses only the function f
# and not its derivative. 
#   
#     dy/dt = f(t,y)
#
# source: https://rosettacode.org/wiki/Runge-Kutta_method
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

""" Sample Function """
def f (t,y):
    return t*np.sqrt(y)
""" Known analitical solution """
def exact_solution (t):
    return (t**2 + 4)**2/16

dt = 0.1
Nstep = 1000
total_time = Nstep*dt
t0 = 0
y0 = 1

def RK4(f, dt, y0, t0=0, Nstep = round(total_time/dt)):
    dy1 = np.zeros(Nstep); dy2 = np.zeros(Nstep); dy3 = np.zeros(Nstep)
    dy4 = np.zeros(Nstep); y = np.zeros(Nstep); t = np.zeros(Nstep)
    t[0] = t0; y[0] = y0
    for i in range(0,Nstep-1):
        dy1[i] = dt*f(t[i], y[i])
        dy2[i] = dt*f(t[i]+0.5*dt, y[i]+0.5*dy1[i])
        dy3[i] = dt*f(t[i]+0.5*dt, y[i]+0.5*dy2[i])
        dy4[i] = dt*f(t[i]+dt, y[i]+dy3[i])
        
        y[i+1] = y[i] + (dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])/6
        t[i+1] = t[i] + dt
    return y 

y = RK4(f, dt, Nstep=Nstep, y0=y0)

time = np.linspace(t0, total_time, round(total_time/dt))
error = np.abs(y - exact_solution(time))

x = np.linspace(0,total_time,1000)
plt.plot(time, y, label='Numerical')
plt.plot(x, exact_solution(x), label='Analitical')
plt.legend()
# %%
""" We analyze the error at the last step
for different values of time step dt
It should behave as error ~ dt^5 """

time_steps = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

error = np.zeros(len(time_steps))
for i in range(0,len(time_steps)):
    y = RK4(f, time_steps[i], Nstep, y0)
    time = np.linspace(t0, total_time, Nstep)
    error[i] = np.abs((exact_solution(time[-1])-exact_solution(time[-2]))-\
               (y[-1]-y[-2]))
               
dt5 = np.power(time_steps,5)               
plt.scatter(dt5, error, label='Local Truncation Error', c='black')

""" Fit polinomiale: verifichiamo pendenza retta in scala log-log """
coef = np.polyfit(dt5, error, 1)
poly = np.poly1d(coef)

xx = np.linspace(min(dt5),max(dt5),1000)
plt.xlim(min(dt5),max(dt5))
plt.xlabel("dt$^5$", fontsize=15)
plt.ylabel("Error", fontsize=15)
plt.plot(xx, poly(xx), 'r--', label='Error $\propto$ dt$^5$') 
plt.legend()

# %%