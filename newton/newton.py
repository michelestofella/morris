# =============================================================================
#
# NEWTON ALGORITHM
#   
# Approximate solution of the equation f(x)=0
#
# =============================================================================

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
    xn: n-th iteration, i.e. the zero of the function
    
    Examples
    --------
    >>> def f(x):
    >>>    return x**2-4
    >>> def Df(x):
    >>>    return 2*x
    >>> x0 = 3.0
    >>> newton(f,Df,x0)
    2.0
    """
    xn = x0
    for n in range(0,max_n):
        f_xn = f(xn)
        if abs(f_xn) < eps:
            # print('Number of iterations to find the solution: ',n)
            return xn
        Df_xn = Df(xn)
        if Df_xn == 0:
            raise Exception('Zero Derivative.\n xn = {}'.format(xn))
            return xn
        xn = xn - f_xn/Df_xn
    raise Exception('Maximum number of iterations reached.\n xn = {}'.format(xn))
    return xn
# %%