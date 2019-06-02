# =============================================================================
#
# RUNGE KUTTA ALGORITHM
# Second Version
#
# This version of the algorithm aims to solve a system of ODEs
#   dy1/dt = f1(t,y1,...,yn)
#   ...
#   dyn/dt = fn(t,y1,...,yn)
#
# =============================================================================

import numpy as np

def f1(t,y1,y2):
    return np.sin(t)+np.cos(y1)+np.sin(y2)
def f2(t,y1,y2):
    return np.cos(t)+np.sin(y2)

""" The following function solves a system of ODEs """
def RK4_system(f, dt, y0, t0, Nstep):
    Neq = len(f)
    t = np.zeros(Nstep+1); t[0] = t0
    y = np.zeros([Neq, Nstep+1]) 
    dy1 = np.zeros([Neq, Nstep+1]); dy2 = np.zeros([Neq, Nstep+1]); 
    dy3 = np.zeros([Neq, Nstep+1]); dy4 = np.zeros([Neq, Nstep+1])
    for j in (0, Neq-1):
        y[j][0] = y0[j]
    
    for i in range(0,Nstep):
        for j in range(0,Neq):
            dy1[j] = dt*f[j](t,*y)
            dy2[j] = dt*f[j](t,*y)
            dy3[j] = dt*f[j](t,*y)
            dy4[j] = dt*f[j](t,*y)
        
            t[i+1] = t[i] + dt
            y[j][i+1] = y[j][i]+(dy1[j][i]+2*dy2[j][i]+2*dy3[j][i]+dy4[j][i])/6
    return t, y

# TEST ONE
# We apply the algorithm to a specific system of functions
# we move forward and backward within the same number of steps
# we expect the final value to be equal to the initial one
# within a tolerance tol
def test_one():
    tol = 1e-6
    def f1(t,y1,y2):
        return np.sin(t)+np.cos(y1)+np.sin(y2)
    def f2(t,y1,y2):
        return np.cos(t)+np.sin(y2)
    fun = [f1,f2]
    dt = 1e-5; y0 = [-1.0,1.0]; t0 = 0.0; Nstep = 1000   
    t_for, y_for = RK4_system(fun, dt=dt, y0=y0, t0=t0, Nstep=Nstep)
    y_last = [y_for[0][-1],y_for[1][-1]]
    t_back, y_back=RK4_system(fun, dt=-dt,y0=y_last,t0=t0+dt*Nstep, Nstep=Nstep)
    assert (np.abs(y_back[0][-1]-y0[0])<tol) &\
           (np.abs(y_back[1][-1]-y0[1])<tol)

# TEST TWO
# We apply the algorithm to a specific system of functions
# we move from a to b, then from b to c and directly from a to c
# we expect the final values to be equal
# within a tolerance tol
def test_two():
    tol = 1e-6    
    def f1(t,y1,y2):
        return np.sin(t)+np.cos(y1)+np.sin(y2)
    def f2(t,y1,y2):
        return np.cos(t)+np.sin(y2)
    fun = [f1,f2]
    dt = 1e-5; y0 = [-1.0,1.0]; t0 = 0.0
    N_ab = 100; N_bc = 50
    t_ab, y_ab = RK4_system(fun, dt, y0, t0, N_ab)
    y_last = [y_ab[0][-1],y_ab[1][-1]]
    t_bc, y_bc = RK4_system(fun, dt, y_last, t0+dt*N_ab, N_bc)
    t_ac, y_ac = RK4_system(fun, dt, y0, t0, N_ab+N_bc)   
    assert (np.abs(y_ac[0][-1]-y_bc[0][-1])<tol) &\
           (np.abs(y_ac[1][-1]-y_bc[1][-1])<tol)   
    
# %%  