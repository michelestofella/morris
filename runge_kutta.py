# =============================================================================
# INTEGRATORE DI RUNGE-KUTTA di ordine 4
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt

# equazione differenziale dy/dt=f(y,t)
def f(y):
    return -y

# To check the derivatives, use symbolic calculus:
# import sympy as sp
# y = sp.symbols('y')
# sp.diff(f(y),y)

def derivata(f,y):
    return -1
def derivata_seconda(f,y):
    return 0
def derivata_terza(f,y):
    return 0

def RK4 (f, y0, dt=0.01, Nstep=1000):
    y = np.zeros(Nstep); y[0]=y0
    for i in range(Nstep-1):
        y[i+1] = y[i] + dt*f(y[i]) + 0.5*dt**2*derivata(f,y[i]) +\
                        dt**3/6*derivata_seconda(f,y[i]) +\
                        dt**4/24*derivata_terza(f,y[i])
    return y

Nstep=1000
dt=0.01
y0=10

x = np.linspace(0,dt*Nstep,Nstep)    
y = RK4(f, y0)    
plt.scatter(x,y)

# %%
# We test the algorithm in the simple case dy/dt=-y
# where the solution is expected to be an exponential decrease
def exact_exp(y0,x):
    return y0*np.exp(-x)
# We test it by changing the integration step delta_t
delta_t = [1., 0.5, 0.1, 0.05, 0.01]
total_time = 5

for dt in delta_t:
    x = np.linspace(0, total_time, round(total_time/dt))
    y = RK4(f,y0,dt,round(total_time/dt))
    plt.plot(x,y,label="dt="+str(dt))

plt.plot(x,exact_exp(y0,x),c='black',label="Exact")
plt.legend()

# %%