# =============================================================================
# RUNGE-KUTTA ALGORITHM
#   
#     dy/dt = f(t,y)
#
# source: https://rosettacode.org/wiki/Runge-Kutta_method
# =============================================================================

import numpy as np

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

# FIRST TEST
# We apply the algorithm to a specific function
# we move forward and backward
# we expect the final value to be equal to the initial one
# within a tolerance tol
def test_one():
    tol = 1e-6
    def f(t,y):
        return t*np.sqrt(y)
    dt = 0.1; y0 = 1.0; t0 = 0.0; Nstep = 100
    t_for, y_for = RK4(f, dt, y0, t0, Nstep)
    y_last = y_for[-1]
    t_back,y_back= RK4(f,-dt, y_last, t0+dt*Nstep, Nstep)
    assert (y_back[-1] < y0+tol) & (y_back[-1] > y0-tol) 
    
# SECOND TEST
# We apply the algorithm to a specific function
# we move from a to b, then from b to c
# and directly from a to c
# we expect the final values to be equal
# within a tolerance tol
def test_two():
    tol = 1e-6
    def f(t,y):
        return t*np.sqrt(y)
    dt = 0.1; y0 = 1.0; t0 = 0.0
    N_ab = 100; N_bc = 50 
    t_ab, y_ab = RK4(f, dt, y0, t0, N_ab)
    t_bc, y_bc = RK4(f, dt, y_ab[-1], t0+N_ab*dt, N_bc)
    t_ac, y_ac = RK4(f, dt, y0, t0, N_ab+N_bc)
    assert (y_ac[-1] < y_bc[-1]+tol) & (y_ac[-1] > y_bc[-1]-tol)

# %%