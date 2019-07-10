import numpy as np
from rk4 import RK4

def test_one():
    def f(t,y):
        return 0
    dt = 0.1; y0 = 1.0; t0 = 0.0; Nstep = 100
    t, y = RK4(f, dt, y0, t0, Nstep)
    assert y[-1]==y0
    
def test_two():
    tol = 1e-6
    def f(t,y):
        return t*np.sqrt(y)
    dt = 0.1; y0 = 1.0; t0 = 0.0; Nstep = 100
    t_for, y_for = RK4(f, dt, y0, t0, Nstep)
    y_last = y_for[-1]
    t_back,y_back= RK4(f,-dt, y_last, t0+dt*Nstep, Nstep)
    assert (y_back[-1] < y0+tol) & (y_back[-1] > y0-tol) 
    
def test_three():
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