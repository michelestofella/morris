# =============================================================================
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

""" The following function solves a system of 2 ODEs """
def RK4_system2(f, dt, y0, t0, Nstep):
    t = np.zeros(Nstep+1); t[0] = t0
    y1 = np.zeros(Nstep+1); y1[0] = y0[0]
    y2 = np.zeros(Nstep+1); y2[0] = y0[1]
    f = [f1,f2]; fun1 = f[0]; fun2 = f[1]

    for i in range(0,Nstep):
        dy1 = [dt*fun1(t,y1,y2), 
               dt*fun2(t,y1,y2)]
        dy2 = [dt*fun1(t+0.5*dt, y1+0.5*dy1[0], y2+0.5*dy1[1]),
               dt*fun2(t+0.5*dt, y1+0.5*dy1[0], y2+0.5*dy1[1])]
        dy3 = [dt*fun1(t+0.5*dt, y1+0.5*dy2[0], y2+0.5*dy2[1]),
               dt*fun2(t+0.5*dt, y1+0.5*dy2[0], y2+0.5*dy2[1])]    
        dy4 = [dt*fun1(t+dt, y1+dy3[0], y2+dy3[1]),
               dt*fun2(t+dt, y1+dy3[0], y2+dy3[1])]    

        t[i+1] = t[i] + dt
        y1[i+1] = y1[i]+(dy1[0][i]+2*dy2[0][i]+2*dy3[0][i]+dy4[0][i])/6
        y2[i+1] = y2[i]+(dy1[1][i]+2*dy2[1][i]+2*dy3[1][i]+dy4[1][i])/6
    y = [y1, y2]
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
    dt = 0.01; y0 = [-1.0,1.0]; t0 = 0.0; Nstep = 1000    
    t_for, y_for = RK4_system2(fun, dt=dt, y0=y0, t0=t0, Nstep=Nstep)
    y_last = [y_for[0][-1],y_for[1][-1]]
    t_back, y_back=RK4_system2(fun, dt=-dt,y0=y_last,t0=t0+dt*Nstep, Nstep=Nstep)
    assert (np.abs(y_back[0][-1]-y0[0])<tol) & (np.abs(y_back[1][-1]-y0[1])<tol)

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
    dt = 0.01; y0 = [-1.0,1.0]; t0 = 0.0
    N_ab = 100; N_bc = 50
    t_ab, y_ab = RK4_system2(fun, dt, y0, t0, N_ab)
    y_last = [y_ab[0][-1],y_ab[1][-1]]
    t_bc, y_bc = RK4_system2(fun, dt, y_last, t0+dt*N_ab, N_bc)
    t_ac, y_ac = RK4_system2(fun, dt, y0, t0, N_ab+N_bc)   
    assert (np.abs(y_ac[0][-1]-y_bc[0][-1])<tol) &\
           (np.abs(y_ac[1][-1]-y_bc[1][-1])<tol)   
    
# %%    