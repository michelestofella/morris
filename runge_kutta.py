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

def RK4(f, dt, y0, Nstep, t0):
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

# %%
    
""" Sample Function 1 """
def f (t,y):
    return t*np.sqrt(y)
""" Solution Function 1 """
def exact_solution (t):
    return (t**2 + 4)**2/16

dt = 0.1
Nstep = 10
t0 = 0.0
y0 = 1.0

y = RK4(f, dt, Nstep=Nstep, y0=y0, t0=t0)

time = np.linspace(t0, Nstep*dt, Nstep)
x = np.linspace(t0, Nstep*dt, 1000)

plt.plot(time, y, label='Numerical')
plt.plot(x, exact_solution(x), label='Analitical')
plt.legend()

# %%

""" We analyze the error at the last step
for different values of time step dt
It should behave as error ~ dt^4 """

time_steps = np.linspace(0.01,1.0,20)

""" We calculate the Global Truncation Error for several dt """
error = np.zeros(len(time_steps))
for i in range(0,len(time_steps)):
    y = RK4(f, time_steps[i], Nstep=Nstep, y0=y0, t0=t0)
    time = np.linspace(t0, Nstep*time_steps[i], Nstep)
    error[i] = np.abs(y[-1]-exact_solution(time[-1]))
                              
dt4 = np.power(time_steps,4)               
plt.scatter(dt4,error)

""" Linear Fit Error vs dt5"""
coef = np.polyfit(dt4, error, 1)
poly = np.poly1d(coef)

xx = np.linspace(min(dt4),max(dt4),1000)
""" We visualize the error with respect to log10(dt^4)
so that we can analyze several order of magnitude for dt """ 
plt.scatter(dt4, error, label='Global Truncation Error', c='black')
plt.plot(xx, poly(xx), 'r--', label='Error $\propto$ $\Delta$t$^4$') 
plt.xlabel("$\Delta$t$^4$", fontsize=15)
plt.ylabel("Error", fontsize=15)
plt.legend()

# %%

""" Sample Function 2 """
def f2(y,t):
    return np.exp(t-y)
""" Solution Function 2 """
def sol_f2(t):
    return t

dt = 0.01
Nstep = 1000

y2 = RK4(f2, dt, Nstep=Nstep, y0=0.0, t0=0.0)

time = np.linspace(t0, t0+Nstep*dt, Nstep)
x = np.linspace(t0, t0+Nstep*dt, 1000)

plt.plot(time, y2, label='Numerical')
plt.plot(x, sol_f2(x), label='Analitical')
plt.legend()

# %%
