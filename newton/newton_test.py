import numpy as np    
from newton import newton
from hypothesis import given
import hypothesis.strategies as st

@given(st.floats(-10,10))
def test_parabola_without_constant_term(s):
    '''
    tests if the Newton algorithm correctly finds that the function 
    f(x)=x^2 has a unique zero at x=0 starting from different initial conditions 
    in between the interval [-10,10].
    '''
    def f(x):
        return x**2
    def Df(x):
        return 2*x
    assert round(newton(f,Df,x0=s),5) == 0

@given(st.floats(-10,10),st.floats(1,16))
def test_parabola_with_constant_term(s,r):
    '''
    tests if the Newton algorithm correctly finds that the function 
    f(x) = x^2-r has two zeros x=sqrt(r) and x=-sqrt(r), with r>0. 
    The strategy is to test the algorithm starting from different initial points: 
    if the initial guess is positive, the positive zero (x=2) is found; 
    if the initial guess is negative, the negative zero (x=-2) is found. 
    '''
    def f(x):
        return x**2-r
    def Df(x):
        return 2*x
    x0 = s
    if x0 > 0:
        assert round(newton(f,Df,x0),5) == round(np.sqrt(r),5)
    elif x0 < 0:
        assert round(newton(f,Df,x0),5) == round(-np.sqrt(r),5)
    
# %%
    
