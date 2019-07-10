import numpy as np    
from newton import newton

def test_one():
    def f(x):
        return x**2-4
    def Df(x):
        return 2*x
    x0_1 = 3.0; x0_2 = -1.5
    assert newton(f,Df,x0_1) == 2.0
    assert newton(f,Df,x0_2) == -2.0

def test_two():
    def f(x):
        return r*x - x**2
    def Df(x):
        return r - 2*x
    r = 3
    assert newton(f,Df,x0=2.0) == 3

def test_three():
    def h(x):
        return np.sin(x)
    def Dh(x):
        return np.cos(x)
    x0 = 0.5
    assert newton(h,Dh,x0) == 0

# %%
    
