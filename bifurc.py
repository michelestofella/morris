# =============================================================================
# 
# EXAMPLES OF BIFURCATION ANALYSIS
#
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt    

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
    start_points = np.linspace(xmin,xmax,100) 
    for x0 in start_points:
        zero.append(round(newton(f,Df,x0,eps),4))
    for element in zero:
        if element not in zeros:
            zeros.append(element)
    return zeros

# %%
""" Monodimensional Case 
1) Saddle Node Bifurcation: 
   dx/dt = phi(x,r) = -r+x**2 """
def f(x):
    return -r+x**2
def Df(x):
    return 2*x

points = []
r_values = np.linspace(0.5,5,100)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-50, xmax=50, eps=1e-12)
    points.append(zeros) 

for i in range(0,len(r_values)):
    plt.scatter(r_values[i],points[i][0],color='red')
    plt.scatter(r_values[i],points[i][1],color='blue')
    
# %%
""" Monodimensional Case
2) Transcritic Bifurcation
    dx/dt = phi(x,r) = r*x-x**2 """
def f(x):
    return r*x-x**2
def Df(x):
    return r-2*x

points = []
r_values = np.linspace(-5,5,100)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-50, xmax=50, eps=1e-12)
    points.append(zeros) 

for i in range(0,len(r_values)):
    plt.scatter(r_values[i],points[i][0],color='red')
    plt.scatter(r_values[i],points[i][1],color='blue')

# %%    
""" Monodimensional Case
3) Fork Bifurcation
    dx/dt = phi(x,r) = r*x-x**3 """
def f(x):
    return r*x-x**3
def Df(x):
    return r-3*x**2
""" Determine equilibrium points for different parameters """
stable_points = []; unstable_points = []
r_values = np.linspace(-5,5,100)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-5, xmax=5, eps=1e-12)
    for i in range(0,len(zeros)):
        if Df(zeros[i])<=0:
            stable_points.append([r,zeros[i],'stable'])
        if Df(zeros[i])>0:
            unstable_points.append([r,zeros[i],'unstable'])
""" Visualization of the results """    
for i in range(0,len(stable_points)):
    plt.scatter(stable_points[i][0],stable_points[i][1],color='red')
for j in range(0,len(unstable_points)):
    plt.scatter(unstable_points[j][0],unstable_points[j][1],color='blue') 
plt.xlabel('parameter r', fontsize=15)
plt.ylabel('Equilibrium Point $x^*$', fontsize=15)
       
# %%   

import pandas as pd

col = ['par','eq_point','stability']
df_stable = pd.DataFrame(stable_points,columns=col)
df_unstable = pd.DataFrame(unstable_points,columns=col)
df = pd.concat([df_stable,df_unstable])

plt.scatter(df['par'],df['eq_point'])

# %%