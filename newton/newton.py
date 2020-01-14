# =============================================================================
#
# NEWTON ALGORITHM
#   
# Approximate solution of the equation f(x)=0
#
# =============================================================================

class Result:
    def __init__(self, x, success, message):
        self.x = x
        self.success = success
        self.message = message

def newton(f, Df, x0, eps=1e-14, max_n=100):
    """ Finds the solution of the equation f(x)=0 
    
    Parameters
    ----------
    f: function to be put in the equation
    Df: derivative of the function 
    x0: initial guess
    eps: precision of the returned value, i.e. how much the function is close to 0 (default: 1e-14)
    max_n: maximum number of iterations (default: 100)
    
    Returns
    -------
    res: Result object with three attributes. ''x'' is the value reached at
         last iteration; ''success'' is a boolean flag indicating if the algorithm
         correctly converged; ''message'' is a string containing information on 
         why the algorithm did not converge.
    
    Examples
    --------
    >>> def f(x):
    >>>    return x**2-4
    >>> def Df(x):
    >>>    return 2*x
    >>> x0 = 3.0
    >>> res = newton(f,Df,x0)
    >>> res.x
    2.0
    """    
    xn = x0
    for n in range(0,max_n):
        f_xn = f(xn)
        if abs(f_xn) < eps:
            # print('Number of iterations to find the solution: ',n)
            success = True
            message = 'Success'
            return Result(xn,success,message)
        Df_xn = Df(xn)
        if Df_xn == 0:
            success = False
            message = 'Zero Derivative.'
            return Result(xn,success,message)
        xn = xn - f_xn/Df_xn
    success = False
    message = 'Maximum number of iterations reached.'
    return Result(xn,success,message)

# %%