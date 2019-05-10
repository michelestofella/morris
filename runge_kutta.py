# =============================================================================
# INTEGRATORE DI RUNGE-KUTTA di ordine 4
# =============================================================================

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# equazione differenziale dy/dt=f(y,t)
def f(y):
    return -y
def derivata(f):
    return -1
def derivata_seconda(f):
    return 0
def derivata_terza(f):
    return 0


def RK4 (f, y0, dt=0.01, Nstep=1000):
    y = np.zeros(Nstep); y[0]=y0
    for i in range(Nstep-1):
        y[i+1] = y[i] + dt*f(y[i]) + 0.5*dt**2*derivata(f(y[i])) +\
                        dt**3/6*derivata_seconda(f(y[i])) +\
                        dt**4/24*derivata_terza(f(y[i]))
    return y

Nstep=1000
dt=0.01

x = np.linspace(0,dt*Nstep,Nstep)    
y = RK4(f, y0=10)    
plt.scatter(x,y)



# %%