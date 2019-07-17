from rk4_system import RK4_system
import numpy as np

def test_1():
    def f(t,y):
        return 0
    fun = [f]; dt = 0.01; y0 = [1.0]; t0 = 0.0; Nstep = 1000
    t, y = RK4_system(fun, dt, y0, t0, Nstep)
    assert y[0][-1] == y0

def test_2():
    def f(t,y):
        return 1
    fun = [f]; dt = 0.01; y0 = [0.0]; t0 = 0.0; Nstep = 1000
    t, y = RK4_system(fun, dt, y0, t0, Nstep)
    assert y[0][-1] == t[-1]
    
def test_3():
    def f1(t,y1,y2):
        return 0
    def f2(t,y1,y2):
        return 1
    fun = [f1,f2]; dt = 0.01; y0 = [0.0,0.0]; t0 = 0.0; Nstep = 1000
    t, y = RK4_system(fun, dt, y0, t0, Nstep)
    assert (y[0][-1] == y0[0]) & (y[1][-1] == t[-1])

def test_4():
    def f1(t,y1,y2):
        return -y1
    def f2(t,y1,y2):
        return -y2
    fun = [f1,f2]; dt = 0.01; y0 = [1.0,1.0]; t0 = 0.0; Nstep = 1000
    t, y = RK4_system(fun, dt, y0, t0, Nstep)
    assert y[0][-1] == y[1][-1]

def test_5():
    def f1(t,y1):
        return -2
    fun = [f1]; dt = 0.001; y0 = [1.0]; t0 = 0.0; Nstep = 1000
    t_for, y_for = RK4_system(fun, dt, y0, t0, Nstep)
    t_back, y_back = RK4_system(fun, -dt, y0=[y_for[0][-1]], t0=t0+Nstep*dt, Nstep=Nstep)
    assert y_back[0][-1] == y_for[0][0]


# %%