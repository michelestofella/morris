# =============================================================================
# 
# Runge Kutta Algorithm
#
# 4th order Runge Kutta algorithm to solve a differential equation
# dy/dt = f(t,x)
#
# =============================================================================

import numpy as np

def RK4(f, dt, y0, t0, Nstep):
    """ 4th order Runge Kutta algorithm
    Solves a differential equation of the form dy/dt = f(t,x) 
    
    Parameters
    ----------
    f: funtion to be integrated, left side of the differential equation
    dt: integration time step 
    y0: initial condition
    t0: initial time 
    Nstep: number of steps to be performed
    
    Returns
    -------
    t: float64 array of size Nstep+1 that contains the time steps
    y: float64 array of size Nstep+1 that contains the solution of the differential equation
    
    Examples
    --------
    >>> def f(t,y):
    >>>    return t*np.sqrt(y)
    >>> dt = 0.01; y0 = 1.0; t0 = 0.0; Nstep = 1000
    >>> RK4(f, dt, y0, t0, Nstep)
    (array([ 0.  ,  0.01,  0.02, ...,  9.98,  9.99, 10.  ]),
     array([  1.        ,   1.00005   ,   1.00020001, ..., 670.81518   ,
            673.4037975 , 675.99999999]))
    """
    y = np.zeros(Nstep+1); t = np.zeros(Nstep+1);
    dy1 = np.zeros(Nstep+1); dy2 = np.zeros(Nstep+1)
    dy3 = np.zeros(Nstep+1); dy4 = np.zeros(Nstep+1)
    t[0] = t0; y[0] = y0
    for i in range(0,Nstep):
        dy1[i] = dt*f(t[i], y[i])
        dy2[i] = dt*f(t[i]+0.5*dt, y[i]+0.5*dy1[i])
        dy3[i] = dt*f(t[i]+0.5*dt, y[i]+0.5*dy2[i])
        dy4[i] = dt*f(t[i]+dt,     y[i]+dy3[i])
        
        y[i+1] = y[i] + (dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])/6
        t[i+1] = t[i] + dt
    return t,y 

# %%