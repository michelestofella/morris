import numpy as np
from rk4 import RK4

def test_one():
    def f(t,y):
        return 0
    dt = 0.1; y0 = 1.0; t0 = 0.0; Nstep = 100
    t, y = RK4(f, dt, y0, t0, Nstep)
    assert y[-1]==y0
    
# %%