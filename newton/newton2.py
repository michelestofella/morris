# =============================================================================
#
# BIDIMENSIONAL NEWTON ALGORITHM
#   
# Approximate solution of the equation system f1(x1,x2)=0, f2(x1,x2)=0
#
# source: http://mathfaculty.fullerton.edu/mathews/n2003/FixPointNewtonMod.html
#
# =============================================================================

import numpy as np
from numpy.linalg import inv
from numpy.linalg import det

class Result:
    def __init__(self, x, success, message):
        self.x = x
        self.success = success
        self.message = message

def newton2(f,Jf,p0,eps=1e-8,max_iter=20):
    """ Bidimensional Newton algorithm: finds the solutions of the system of equations
    f1(x,y) = 0
    f2(x,y) = 0
    
    Parameters
    ----------
    f: list of the two functions that form the system of equations
    Jf: jacobian matrix of f
    p0: initial guess, bidimensional list
    eps: precision of the returned value, i.e. how close f is close to zero (default: 1e-8)
    max_iter: maximum number of iterations
    
    Returns
    -------
    res: Result object with 3 attributes: ''x'' is a bidimensional array containing
       the zero of the function, ''success'' is a boolean flag indicating if 
       the algorithm correctly coverged; ''message'' is a string containing information
       on why the algorithm did not converge properly. 
        
    Examples
    --------
    >>> p0 = [0.1, 0.7]
    >>> def f1(x,y):
    >>>     return 1 - 4*x + 2*x**2 - 2*y**3
    >>> def f2(x,y):
    >>>     return -4 + x**4 + 4*y + 4*y**4
    >>> def f(x,y):
    >>>     return [f1(x,y),f2(x,y)]
    >>> def Jf(x,y):
    >>>     return [[-4+4*x, -6*y**2],
    >>>             [4*x**3, 4+16*y**3]]
    >>> res = newton2(f=f,Jf=Jf,p0=p0)
    >>> res.x
    array([0.06177013, 0.72449052])
    """
    for k in range(0,max_iter):
        f_k = f(p0[0],p0[1])
        Jf_k = Jf(p0[0],p0[1])
        if det(Jf_k) == 0:
            if k == 0:
                pk = p0
            success = False
            message = 'Zero determinant'
            return Result(pk,success,message)

        invJf_k = inv(Jf_k)
        pk = p0 - np.dot(invJf_k,f_k)
    
        dist = np.sqrt(np.dot(pk-p0,pk-p0))
        if dist < eps:
            success = True
            message = 'Success'
            return Result(pk,success,message)
        p0 = pk
    success = False
    message = 'Max number of iterations reached'
    return Result(pk,success,message)
    
# %%