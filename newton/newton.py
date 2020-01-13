# =============================================================================
#
# NEWTON ALGORITHM
#   
# Approximate solution of the equation f(x)=0
#
# =============================================================================

def newton(f, Df, x0, eps=1e-16, max_n=50):
    """ Finds the solution of the equation f(x)=0 
    
    Parameters
    ----------
    f: function to be put in the equation
    Df: derivative of the function 
    x0: initial guess
    eps: precision of the returned value, i.e. how much the function is close to 0 (default: 1e-16)
    max_n: maximum number of iterations (default: 50)
    
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
    Number of iterations to find the solution: 5
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
            print('Zero derivative')
            return None
        xn = xn - f_xn/Df_xn
    print('Maximum number of iterations reached')
    return xn
 
# %%