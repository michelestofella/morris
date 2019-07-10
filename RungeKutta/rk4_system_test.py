from rk4_system import RK4_system
import numpy as np

""" Tests for one ODE """

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
 
# %%
    
""" Test for a system of two ODEs """

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

# %%

""" More -cryptic- tests """
def test_one():
    tol = 1e-6
    def f(t,y):
        return t*np.sqrt(y)
    dt = 0.1; y0 = [1.0]; t0 = 0.0; Nstep = 100
    fun = [f]
    t_for, y_for = RK4_system(fun, dt, y0, t0, Nstep)
    y_last = [y_for[0][-1]]
    t_back,y_back= RK4_system(fun,-dt, y_last, t0+dt*Nstep, Nstep)
    assert (np.abs(y_back[0][-1]-y0[0] < tol))

def test_two():
    tol = 1e-6
    def f(t,y):
        return t*np.sqrt(y)
    dt = 0.1; y0 = [1.0]; t0 = 0.0
    N_ab = 100; N_bc = 50 
    fun = [f]
    t_ab, y_ab = RK4_system(fun, dt, y0, t0, N_ab)
    t_bc, y_bc = RK4_system(fun, dt, [y_ab[0][-1]], t0+N_ab*dt, N_bc)
    t_ac, y_ac = RK4_system(fun, dt, y0, t0, N_ab+N_bc)
    assert (np.abs(y_ac[0][-1]-y_bc[0][-1] < tol))

# %%

def test_three():
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

def test_four():
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

def test_five():
    tol = 1e-6
    def fun1(t,x,y,z):
        return -x+3*z
    def fun2(t,x,y,z):
        return -y+2*z
    def fun3(t,x,y,z):
        return x**2-2*z
    f = [fun1, fun2, fun3]
    dt = 1e-6; y0 = [0.0, .5, 3.0]; t0 = 0.0; Nstep = 1000

    t_for, y_for = RK4_system(f=f, dt=dt, y0=y0, t0=t0, Nstep=Nstep)
    y_last = [y_for[0][-1], y_for[1][-1], y_for[2][-1]]
    t_back,y_back= RK4_system(f=f,dt=-dt,y0=y_last,t0=t0+dt*Nstep,Nstep=Nstep)
    assert (np.abs(y_back[0][-1]-y0[0])<tol) &\
           (np.abs(y_back[1][-1]-y0[1])<tol) &\
           (np.abs(y_back[2][-1]-y0[2])<tol)

# %%