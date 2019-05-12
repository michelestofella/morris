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
        y[i+1] = y[i] + dt*f(y[i]+0.5*dt*f(y[i])) +\
                        dt**2/2*derivata(f,y[i]) +\
                        dt**3/6*derivata_seconda(f,y[i]) +\
                        dt**4/24*derivata_terza(f,y[i])
    return y

Nstep=1000
dt=0.01
y0=10

x = np.linspace(0,dt*Nstep,Nstep)    
y = RK4(f, y0)    
plt.plot(x,y)

# %%
# We test the algorithm in the simple case dy/dt=-y
# where the solution is expected to be an exponential decrease
def exact_exp(x,y0):
    return y0*np.exp(-x)
# We test it by changing the integration step delta_t
delta_t = [0.1, 0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]
total_time = 5

for dt in delta_t:
    x = np.linspace(0, total_time, round(total_time/dt))
    y = RK4(f,y0,dt,round(total_time/dt))
    plt.plot(x,y,label="dt="+str(dt))

plt.plot(np.linspace(0,5,1000),exact_exp(np.linspace(0,5,1000),y0),c='black',label="Exact")
plt.legend()

# %%
# Second part of the test: we evaluate the distance between the 
# approximated curves and the exact curve

# Error between the dt=1.00 curve and the exact one
# x = np.linspace(0, total_time, round(total_time/delta_t[0]))
# y = RK4(f,y0,delta_t[0],round(total_time/delta_t[0]))
# exact_y = exact_exp(x,y0)

# err = np.zeros(round(total_time/delta_t[0]))
# for i in range(0,round(total_time/delta_t[0])):
#     err[i] = np.abs(y[i]-exact_y[i])
# plt.scatter(x,err)

# I want a plot Error vs Deltat
mean_err = np.zeros(len(delta_t))
for i in range(0,len(delta_t)):
    x = np.linspace(0, total_time, round(total_time/delta_t[i]))
    y = RK4(f,y0,delta_t[i],round(total_time/delta_t[i]))
    exact_y = exact_exp(x,y0)
    mean_err[i] = np.abs(y[round(total_time/delta_t[i])-1]-exact_y[round(total_time/delta_t[i])-1])

plt.scatter(delta_t,mean_err, c='black')
plt.xlim(0,1.0); plt.ylim(0,1.0)
plt.xlabel("dt", fontsize=18); plt.ylabel("Mean Error", fontsize=18)

coef = np.polyfit(delta_t, mean_err, 2)
poly = np.poly1d(coef)
xx = np.linspace(0,1.0,1000)
plt.plot(xx, poly(xx), 'r--') 
   
# %%
# =============================================================================
# SECOND VERIONS OF RUNGE-KUTTA ALGORITHM
# =============================================================================

# http://www2.hawaii.edu/~jmcfatri/math407/RungeKuttaTest.html










# %%