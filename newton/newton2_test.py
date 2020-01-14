from newton2 import newton2
import numpy as np
import pytest
from hypothesis import given
import hypothesis.strategies as st

@given(st.floats(-10,10),st.floats(-10,10))
def test_unique_solution(x0,y0):
    ''' 
    tests if the algorithm correctly finds that the unique zero of the system
      f(x,y) = x = 0
      g(x,y) = y = 0
    is (x,y)=(0,0). The strategy is to use different initial guesses for both x and y.
    '''
    def f1(x,y):
        return x
    def f2(x,y):
        return y
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[1,0],
                [0,1]]
    p0 = [x0,y0]
    xn, yn = newton2(f,Jf,p0).x
    assert round(xn,6) == 0
    assert round(yn,6) == 0

@given(st.floats(-10,10),st.integers(-10,10),st.integers(1,10))
def test_two_possible_solutions(x0,y0,r):
    '''
    tests if the algorithm finds that the zeros of the system
      f1(x,y) = x-y = 0
      f2(x,y) = y**2-r = 0
    has two solutions at (x,y)=(-sqrt(r),-sqrt(r)) and (x,y)=(sqrt(r),sqrt(r)).
    The strategy is to test the outcome starting from different initial guesses
    and using different values of the parameter r in the function f2(x,y).
    Notice that if the initial guess on y is y0=0, the determinant of the jacobian 
    is zero, thus the algorithm should raise an exception.
    '''
    def f1(x,y):
        return x - y
    def f2(x,y):
        return y**2 - r
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[1,1],
                [0,2*y]]
    p0 = [x0,y0]
    ''' if y0=0, the determinant of the jacobian is zero,
    thus an exception should be raised. '''
    if y0 == 0:
        assert newton2(f,Jf,p0).success == False
    else:
        xn, yn = newton2(f,Jf,p0).x
        assert round(xn,6) == round(yn,6)
        if y0 > 0:
            assert round(yn,6) == round(np.sqrt(r),6)
        if y0 < 0:
            assert round(yn,6) == round(-np.sqrt(r),6)        
        
@given(st.floats(-10,10),st.floats(-10,10),
       st.integers(-10,10),st.integers(-10,10))
def test_zero_determinant_exception(x0,y0,
               r1,r2):
    '''
    tests if in a particular case where the determinant of the jacobian is
    zero, namely the system
      f1(x,y) = r1
      f2(x,y) = r2
    where r1,r2 are constants, the algorithm raises an exception. 
    The strategy is to check the algorithm at different constant values r1,r2 
    and starting from different initial guesses.
    '''
    p0 = [x0,y0]
    def f1(x,y):
        return r1
    def f2(x,y):
        return r2
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[0,0],
                [0,0]]
    assert newton2(f,Jf,p0).success == False
    
@given(st.floats(-10,10),st.floats(-10,10))
def test_max_iterations_exception(x0,y0):
    '''
    tests if the maximum number of iterations is reached in a system
      f1(x,y) = x^2 + y^2 + 1
      f2(x,y) = 2x
    where no zero exist. The strategy is to check the algorithm starting
    from different initial guesses. 
    '''
    def f1(x,y):
        return x**2+y**2+1
    def f2(x,y):
        return 2*x
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[2*x,2*y],
                [2,0]]
    p0 = [x0,y0]
    assert newton2(f,Jf,p0).success == False

# %%