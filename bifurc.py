# =============================================================================
# 
# EXAMPLES OF BIFURCATION ANALYSIS
#
# =============================================================================

""" Newton Algorithm: to find a solution to the equation f(x)=0 """
def newton(f, Df, x0, eps, max_n=50):
    xn = x0
    for n in range(0,max_n):
        f_xn = f(xn)
        if abs(f_xn) < eps:
            return xn
        Df_xn = Df(xn)
        if Df_xn == 0:
            print('Zero derivative')
            return None
        xn = xn - f_xn/Df_xn
    print('Maximum number of iterations reached')
    return None

""" Newton algorithm applied to a set of initial points: from xmin to xmax """
def find_zeros(f, Df, xmin, xmax, eps, max_n=50):
    zero = []; zeros = []
    start_points = np.linspace(xmin,xmax,20) 
    for x0 in start_points:
        zero.append(round(newton(f,Df,x0,eps),4))
    for element in zero:
        if element not in zeros:
            zeros.append(element)
    return zeros

# %%

import numpy as np
import matplotlib.pyplot as plt    

""" Monodimensional Case 
1) Saddle Node Bifurcation: dx/dt = phi(x,r) = -r+x**2 """
def f(x):
    return -r+x**2
def Df(x):
    return 2*x

points = []
r_values = np.linspace(0.5,5,10)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-50, xmax=50, eps=1e-5)
    points.append(zeros) 

for i in range(0,len(r_values)):
    plt.scatter(r_values[i],points[i][0],color='red')
    plt.scatter(r_values[i],points[i][1],color='blue')
    
# %%