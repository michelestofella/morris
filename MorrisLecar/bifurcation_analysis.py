# =============================================================================
# 
# SCRIPT FOR BIFURCATION ANALYSIS OF THE MORRIS LECAR MODEL
# Bifurcation analysis is performed
# Bifurcation diagram is generated
# 
# Several parameters must be parsed in order to run the algorithm:
#   --v_ca  the parameter of the morris lecar model that discriminates
#           different classes of neurons.
#   --Imin  minimum value of applied current analysed
#   --Imax  maximum value of applied current analysed 
#   --v0min minimum value of initial condition on voltage
#   --v0max maximum value of initial condition on voltage
#   --out   name of the generated figure
#
# =============================================================================

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--v_ca")
parser.add_argument("--Imin")
parser.add_argument("--Imax")
parser.add_argument("--v0min")
parser.add_argument("--v0max")
parser.add_argument("--out")

config = {}
opts = parser.parse_args()

if opts.v_ca:
    v_ca = float(opts.v_ca)
if opts.Imin:
    Imin = float(opts.Imin)
if opts.Imax:
    Imax = float(opts.Imax)
if opts.v0min:
    v0min = float(opts.v0min)
if opts.v0max:
    v0max = float(opts.v0max)
    
if opts.out:
    out = opts.out

# %%

from fixed_parameters import *

""" Model for Bifurcation Analysis 
Here we need to remove the temporal dependency """
def f1(x,y):
    part1 = -g_ca*m_inf(x)*(x-E_ca)
    part2 = -g_k*y*(x-E_k)
    part3 = -g_leak*(x-E_leak)
    return (part1+part2+part3+I_app)/c
def f2(x,y):
    return phi_w*((w_inf(x)-y)/tau_w(x))    

def df1dx(x,y):
    part1 = -g_ca*(x-E_ca)/(2*theta_ca*(np.cosh((v_ca-x)/theta_ca))**2)
    part2 = -g_ca*m_inf(x)
    part3 = -g_k*y
    part4 = -g_leak
    return (part1+part2+part3+part4)/c
def df1dy(x,y):
    return -g_k*(x-E_k)/c

def m_inf(v):
    return 0.5*(1+np.tanh((v-v_ca)/theta_ca))
def w_inf(v):
    return 0.5*(1+np.tanh((v-v_k)/theta_k))
def tau_w(v):
    return 1/(np.cosh((v-v_k)/(2*theta_k)))

def df2dx(x,y):
    part1 = 1/(2*theta_k*(np.cosh((v_k-x)/theta_k))**2)
    part2 = np.sinh((v_k-x)/(2*theta_k))/(theta_k*(np.cosh((v_k-x)/theta_k)+1))
    return phi_w*(part1*tau_w(x)-(w_inf(x)-y)*part2)/(tau_w(x))**2
def df2dy(x,y):
    return -phi_w/tau_w(x)

def f(x,y):
    return [f1(x,y),f2(x,y)]
def Jf(x,y):
    return [[df1dx(x,y),df1dy(x,y)],
            [df2dx(x,y),df2dy(x,y)]]

# %%

""" Set I_app values and v0_values to draw bifurcation diagram """
I_app_values = np.linspace(Imin,Imax,200)    
v0_values = np.linspace(v0min,v0max,61)
w0 = 0.0

""" Find the zeros of the function """
v_zeros = []; w_zeros = []
for I_app in I_app_values:
    zero_v = []; zero_w = []
    zeros_v = []; zeros_w = []
    
    for v0 in v0_values:
        p0 = [v0,w0]
        res = newton2(f,Jf,p0)
        if res.success == True:
            point = res.x
            zero_v.append(round(point[0],5))
            zero_w.append(round(point[1],5))
            for element in zero_v:
                if element not in zeros_v:
                    zeros_v.append(element)
            for element in zero_w:
                if element not in zeros_w:
                    zeros_w.append(element)

    v_zeros.append(zeros_v); w_zeros.append(zeros_w)

# %%

""" Discriminate stability of each zero found previously """
stable = []; unstable = []
for i in range(0,len(I_app_values)):
    I_app = I_app_values[i]
    for j in range(0,len(v_zeros[i])):
        eigen = np.linalg.eigvals(Jf(v_zeros[i][j],w_zeros[i][j]))
        for k in range(0,len(eigen)):
            if eigen[k]>0:
                unstable.append([I_app_values[i],v_zeros[i][j],w_zeros[i][j],'unstable'])
                break
            stable.append([I_app_values[i],v_zeros[i][j],w_zeros[i][j],'stable'])
            break
            
# %%
 
""" Reorganize data """
col_st = ['I_app','V','w','Stability']
stab = pd.DataFrame(stable,columns=col_st)
unstab = pd.DataFrame(unstable, columns=col_st)

unstab_up = unstab.loc[unstab['V']>-30]
unstab_down = unstab.loc[unstab['V']<-30]

# %%

""" Bifurcation Diagram """
plt.figure(figsize=(15,10))
plt.plot(stab['I_app'],stab['V'],'b')
plt.plot(unstab_up['I_app'],unstab_up['V'],'b--')
plt.plot(unstab_down['I_app'],unstab_down['V'],'b--')
plt.ylim(v0min,v0max)
plt.xlabel('$I_{app}$',fontsize=18)
plt.ylabel('V',fontsize=18)
plt.grid(linestyle=':')
plt.savefig(out+'.png')

# %%