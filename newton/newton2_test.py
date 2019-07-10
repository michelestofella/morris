from newton2 import newton2

def test_one():
    p0 = [0.1, 0.7]
    def f1(x,y):
        return 1 - 4*x + 2*x**2 - 2*y**3
    def f2(x,y):
        return -4 + x**4 + 4*y + 4*y**4
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[-4+4*x, -6*y**2],
                [4*x**3, 4+16*y**3]]
    real_sol = [0.061770,0.724491]
    num_sol = newton2(f=f,Jf=Jf,p0=p0)
    assert round(num_sol[0],6) == real_sol[0]
    assert round(num_sol[1],6) == real_sol[1]
    
def test_two():
    p0 = [0.5, -2.0]
    def f1(x,y):
        return x
    def f2(x,y):
        return y
    def f(x,y):
        return [f1(x,y),f2(x,y)]
    def Jf(x,y):
        return [[1, 0],
                [0, 1]]
    num_sol = newton2(f,Jf,p0)
    assert num_sol[0] == 0.
    assert num_sol[1] == 0.

# %%