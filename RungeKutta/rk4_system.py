# =============================================================================
#
# RUNGE KUTTA ALGORITHM
#
# 4th order Runge Kutta algorithm to solve a system of ODEs
#   dy1/dt = f1(t,y1,...,yn)
#   ...
#   dyn/dt = fn(t,y1,...,yn)
# 
# =============================================================================

import numpy as np

def RK4_system(f, dt, y0, t0, Nstep):
    """ 4th order Runge Kutta algorithm
    Solves a set of ODEs of the form 
    dy1/dt = f1(t,y1,...,yn)
    ...
    dyn/dt = fn(t,y1,...,yn)
    
    Parameters
    ----------
    f: list of functions to be integrated (left side vector of the system of ODEs)
    dt: integration time step 
    y0: list of initial conditions
    t0: initial time 
    Nstep: number of steps to be performed
    
    Returns
    -------
    t: float64 array of size Nstep+1 that contains the time steps
    y: float64 array of size (len(f),Nstep+1) that contains the solutions of the system of ODEs
    
    Examples
    --------
    >>> def f1(t,y1,y2):
    >>>     return 0
    >>> def f2(t,y1,y2):
    >>>     return 1
    >>> f = [f1,f2]
    >>> y0 = [0.,0.]; Nstep = 100; dt = 0.01; t0 = 0.
    >>> t, y = RK4_system(f,dt,y0,t0,Nstep)
    >>> t
    array([0.  , 0.01, 0.02, 0.03, ..., 0.96, 0.97, 0.98, 0.99, 1.  ])
    >>> y[0]
    array([0., 0., 0., 0., 0., 0., ..., 0., 0., 0., 0., 0., 0., 0.])
    >>> y[1]
    array([0.  , 0.01, 0.02, 0.03, ..., 0.96, 0.97, 0.98, 0.99, 1.  ])
    """

    Neq = len(f)
    t = np.zeros(Nstep+1)
    t[0] = t0
    y = np.zeros([Neq, Nstep+1]) 
    dy1 = np.zeros([Neq, Nstep+1]); dy2 = np.zeros([Neq, Nstep+1]); 
    dy3 = np.zeros([Neq, Nstep+1]); dy4 = np.zeros([Neq, Nstep+1])
    for j in range(0, Neq):
        y[j][0] = y0[j]
    
    for i in range(0,Nstep):
        for j in range(0,Neq):
            dy1[j] = dt*f[j](t,*y)
            dy2[j] = dt*f[j](t+0.5*dt,*y+0.5*dy1)
            dy3[j] = dt*f[j](t+0.5*dt,*y+0.5*dy2)
            dy4[j] = dt*f[j](t+dt,*y+dy3)
        
            t[i+1] = t[i] + dt
            y[j][i+1] = y[j][i]+(dy1[j][i]+2*dy2[j][i]+2*dy3[j][i]+dy4[j][i])/6
    return t, y

# %%