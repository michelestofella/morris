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

# Sample function
def f (t,y):
    return t*np.sqrt(y)
# Known analitical solution
def exact_solution (t):
    return (t**2 + 4)**2/16

dt = 0.1
total_time = 10
t0 = 0
y0 = 1

def RK4(f, dt, total_time, y0, t0=0):
    Nstep = round(total_time/dt)
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

y = RK4(f, dt, total_time, y0)

time = np.linspace(t0, total_time, round(total_time/dt))
error = np.abs(y - exact_solution(time))

x = np.linspace(0,10,1000)
plt.plot(time, y, label='Numerical')
plt.plot(x, exact_solution(x), label='Analitical')
plt.legend()
# %%
d = {'Numerical' : y, 
     'Analytical' : exact_solution(time),
     'Error' : error}
res = pd.DataFrame(d)
res.head(10)

coef = np.polyfit(time, error, 5)
poly = np.poly1d(coef)
xx = np.linspace(0,10,1000)

plt.scatter(time,error)
plt.plot(xx, poly(xx), 'r--') 

# %%