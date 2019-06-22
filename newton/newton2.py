# =============================================================================
#
# BIDIMENSIONAL NEWTON ALGORITHM
#   
# Approximate solution of the equation system f1(x1,x2)=0, f2(x1,x2)=0
#
# source: http://mathfaculty.fullerton.edu/mathews/n2003/FixPointNewtonMod.html
# =============================================================================

import numpy as np
from numpy.linalg import inv

def newton2(f,Jf,p0,eps=1e-8,max_iter=20):
    for k in range(0,max_iter):
        f_k = f(p0[0],p0[1])
        Jf_k = Jf(p0[0],p0[1])
        invJf_k = inv(Jf_k)
        pk = p0 - np.dot(invJf_k,f_k)
    
        dist = np.sqrt(np.dot(pk-p0,pk-p0))
        if dist < eps:
            print('Number of iterations: ',k)
            print('Solution: ',pk)
            return pk
        p0 = pk
    print('Maximum number of iterations reached: ',max_iter)
    
# %%

def test_one():
    tol = 1e-6
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
    real_sol = [0.0617701,0.724491]
    num_sol = newton2(f=f,Jf=Jf,p0=p0)
    diff = np.abs(num_sol-real_sol)
    assert np.all(diff < tol)

# %%