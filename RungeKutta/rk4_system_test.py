from rk4_system import RK4_system
from hypothesis import given
import hypothesis.strategies as st

@given(st.floats(-10,10),st.integers(10,1000))
def test_monodimensional_constant_case(y0,Nstep):
    ''' 
    tests if in the monodimensional example
      dy/dt = 0
    the solution is a constant, checking that the final value is equal to 
    the initial condition. The strategy is to set different initial conditions
    and to change the total number of steps.
    '''
    def f(t,y):
        return 0
    fun = [f]; dt = 0.01; p0 = [y0]; t0 = 0.0
    t, y = RK4_system(fun, dt, p0, t0, Nstep)
    assert round(y[0][-1],5) == round(y0,5)

@given(st.floats(-10,10),st.integers(10,1000))
def test_monodimensional_linear_case(y0,Nstep):
    '''
    tests if in the monodimensional example
      dy/dt = 1
    the solution is a line y = t + y0. The strategy is to check the algorithm 
    with different initial conditions and with different number of steps.
    '''
    def f(t,y):
        return 1
    fun = [f]; dt = 0.01; p0 = [y0]; t0 = 0.0
    t, y = RK4_system(fun, dt, p0, t0, Nstep)
    assert round(y[0][-1],5) == round(y0+t[-1],5)

@given(st.floats(-10,10),st.floats(-10,10),st.integers(10,1000))
def test_bidimensional_independent_case(x0,y0,Nstep):
    '''
    tests if in the bidimensional independent case
      dx/dt = 0
      dy/dt = 1
    the solution for x is a constant equal to the initial condition x0 and 
    the solution for y is a line y = t + y0. The strategy is to test both the
    asserts at different initial conditions for both x0 and y0, integrating
    runge kutta for different number of steps. 
    '''
    def f1(t,y1,y2):
        return 0
    def f2(t,y1,y2):
        return 1
    fun = [f1,f2]; dt = 0.01; p0 = [x0,y0]; t0 = 0.0
    t, y = RK4_system(fun, dt, p0, t0, Nstep)
    assert round(y[0][-1],5) == round(x0,5)
    assert round(y[1][-1],5) == round(y0+t[-1],5)

@given(st.floats(-10,10),st.integers(10,1000))
def test_bidimensional_dependent_case(x0,Nstep):
    '''
    tests if in the bidimensional case
      dx/dt = -x
      dy/dt = -x
    the solution for x and y is the same (the final value is checked). The 
    strategy is to set different initial conditions on x0 and y0, being sure that 
    x0 = y0 and to integrate the system for a different number of steps.
    '''
    def f1(t,y1,y2):
        return -y1
    def f2(t,y1,y2):
        return -y1
    fun = [f1,f2]; dt = 0.01; y0 = [x0,x0]; t0 = 0.0
    t, y = RK4_system(fun, dt, y0, t0, Nstep)
    assert round(y[0][-1],5) == round(y[1][-1],5)

@given(st.floats(-10,10),st.integers(10,1000),st.integers(1,10))
def test_reversibility(y0,Nstep,r):
    '''
    tests the reversibility of the algorithm in the linear case
      dy/dt = -r
    where r is a constant. 
    1) The algorithm is applied for Nstep with time step dt; the solution goes
    from y0 to yn;
    2) the algorithm is applied backward again for Nstep with time step -dt; the 
    solution goes from yn to y0.
    The test checks if the final value of the backward algorithm is the initial 
    value of the forward algorithm. 
    The strategy is to use different initial conditions (on the forward algorithm),
    integrating the system with different number of steps and varying the 
    constant r of the system. 
    '''
    def f1(t,y1):
        return -r
    fun = [f1]; dt = 0.001; p0 = [y0]; t0 = 0.0
    t_for, y_for = RK4_system(fun, dt, p0, t0, Nstep)
    last_time = t0+Nstep*dt
    last_y = y_for[0][-1]
    t_back, y_back = RK4_system(fun,-dt,[last_y],last_time,Nstep)
    assert round(y_back[0][-1],5) == round(y_for[0][0],5)

# %%