# =============================================================================
# 
# EXAMPLES OF MONODIMENSIONAL BIFURCATION ANALYSIS
#
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt 
import pandas as pd   

import sys
sys.path.insert(0, '../newton')
from newton import newton

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

""" Saddle Node Bifurcation: dx/dt = phi(x,r) = -r+x**2 """
def f(x):
    return -r+x**2
def Df(x):
    return 2*x
        
""" We determine equilibrium points (stable and unstable) 
    for different values of the parameter r """    
stable_points = []; unstable_points = []
r_values = np.linspace(0.1,5,100)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-5, xmax=5, eps=1e-12)
    for i in range(0,len(zeros)):
        if Df(zeros[i])<=0:
            stable_points.append([r,zeros[i],'stable'])
        if Df(zeros[i])>0:
            unstable_points.append([r,zeros[i],'unstable'])  

""" Data organization in a dataframe """
col = ['par','eq_point','stability']
df_stable = pd.DataFrame(stable_points,columns=col)
df_unstable = pd.DataFrame(unstable_points,columns=col)
df = pd.concat([df_stable,df_unstable])

""" Data Visualization """
plt.plot(df_stable['par'],df_stable['eq_point'],
            color='black',label='stable')
plt.plot(df_unstable['par'],df_unstable['eq_point'],
            '--',color='black',label='unstable')
plt.xlabel('Parameter r',fontsize=15)
plt.ylabel('Equilibrium Point $x^*$',fontsize=15)
plt.legend()
 
# %%

""" Transcritic Bifurcation: dx/dt = phi(x,r) = r*x-x**2 """
def f(x):
    return r*x-x**2
def Df(x):
    return r-2*x

""" We determine equilibrium points (stable and unstable) 
    for different values of the parameter r """    
stable_points = []; unstable_points = []
r_values = np.linspace(-5,5,100)
for r in r_values:
    zeros = find_zeros(f, Df, xmin=-5, xmax=5, eps=1e-12)
    for i in range(0,len(zeros)):
        if Df(zeros[i])<=0:
            stable_points.append([r,zeros[i],'stable'])
        if Df(zeros[i])>0:
            unstable_points.append([r,zeros[i],'unstable'])  

""" Data organization in a dataframe """
col = ['par','eq_point','stability']
df_stable = pd.DataFrame(stable_points,columns=col)
df_unstable = pd.DataFrame(unstable_points,columns=col)
df = pd.concat([df_stable,df_unstable])

""" Data Visualization """
plt.plot(df_stable['par'],df_stable['eq_point'],
            color='black',label='stable')
plt.plot(df_unstable['par'],df_unstable['eq_point'],
            '--',color='black',label='unstable')
plt.xlabel('Parameter r',fontsize=15)
plt.ylabel('Equilibrium Point $x^*$',fontsize=15)
plt.legend()

# %%

""" Fork Bifurcation: dx/dt = phi(x,r) = r*x-x**3 """
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

""" Organization of the results in a dataframe """    
col = ['par','eq_point','stability']
df_stable = pd.DataFrame(stable_points,columns=col)
df_unstable = pd.DataFrame(unstable_points,columns=col)
df = pd.concat([df_stable,df_unstable])

""" Data Visualization """
plt.scatter(df_stable['par'],df_stable['eq_point'],s=5,
            color='red',label='stable')
plt.scatter(df_unstable['par'],df_unstable['eq_point'], s=5,
            color='blue',label='unstable')
plt.xlabel('Parameter r',fontsize=15)
plt.ylabel('Equilibrium Point $x^*$',fontsize=15)
plt.legend()

# %%