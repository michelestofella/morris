# =============================================================================
# 
# SCRIPT FOR BIFURCATION ANALYSIS OF THE MORRIS LECAR MODEL
# Bifurcation analysis is performed
# Bifurcation diagram is generated
# 
# =============================================================================

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--v_ca")
parser.add_argument("--dt")
parser.add_argument("--Nstep")
parser.add_argument("--v0")
parser.add_argument("--w0")

parser.add_argument("--Imin")
parser.add_argument("--Imax")
parser.add_argument("--v0min")
parser.add_argument("--v0max")

config = {}
opts = parser.parse_args()

if opts.v_ca:
    v_ca = float(opts.v_ca)
if opts.dt:
    dt = float(opts.dt)
if opts.Nstep:
    Nstep = int(opts.Nstep)
if opts.v0:
    v0 = float(opts.v0)
if opts.w0:
    w0 = float(opts.w0)

if opts.Imin:
    Imin = float(opts.Imin)
if opts.Imax:
    Imax = float(opts.Imax)
if opts.v0min:
    v0min = float(opts.v0min)
if opts.v0max:
    v0max = float(opts.v0max)

# %%

from fixed_parameters import *

""" Model for the Runge Kutta Integration 
Here we need to insert the temporal dependency """
def g1(t,y1,y2):
    part1 = g_ca*m_inf(y1)*(E_ca-y1)
    part2 = g_k*y2*(E_k-y1)
    part3 = g_leak*(E_leak-y1)
    return (part1+part2+part3+I_app)/c
def g2(t,y1,y2):
    return phi_w*((w_inf(y1)-y2)/tau_w(y1))

def m_inf(v):
    return 0.5*(1+np.tanh((v-v_ca)/theta_ca))
def w_inf(v):
    return 0.5*(1+np.tanh((v-v_k)/theta_k))
def tau_w(v):
    return 1/(np.cosh((v-v_k)/(2*theta_k)))

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

count_warn = 0
""" Find the zeros of the function """
v_zeros = []; w_zeros = []
for I_app in I_app_values:
    zero_v = []; zero_w = []
    zeros_v = []; zeros_w = []
    
    for v0 in v0_values:
        p0 = [v0,w0]
        try:
            point = newton2(f,Jf,p0)
            zero_v.append(round(point[0],5))
            zero_w.append(round(point[1],5))
            for element in zero_v:
                if element not in zeros_v:
                    zeros_v.append(element)
            for element in zero_w:
                if element not in zeros_w:
                    zeros_w.append(element)

        except:
            count_warn += 1
    
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

""" Integrate the model at different values of I_app 
    to find maximum and minimum values of the limit cycle """
g = [g1,g2]
t0 = 0.0
y0 = [v0,w0]

max_v = []; min_v = []; freq = []
for j in range(0,len(I_app_values)):
    #print("Calculating:",round(j*100/len(I_app_values),1),"%")
    I_app = I_app_values[j]
    time, sol = RK4_system(g, dt, y0, t0, Nstep)
    if I_app > stable[-1][0]:
        max_v.append([I_app,max(sol[0])])
        min_v.append([I_app,min(sol[0])])
#print('Done: 100.0 %')
            
# %%            
 
""" Reorganize data """
col_st = ['I_app','V','w','Stability']
stab = pd.DataFrame(stable,columns=col_st)
unstab = pd.DataFrame(unstable, columns=col_st)

unstab_up = unstab.loc[unstab['V']>-30]
unstab_down = unstab.loc[unstab['V']<-30]

col=['I_app','V']
vmax = pd.DataFrame(max_v,columns=col)
vmin = pd.DataFrame(min_v,columns=col)

# %%

""" Bifurcation Diagram """
plt.plot(stab['I_app'],stab['V'],'b')
#plt.plot(unstab['I_app'],unstab['V'],'b--')
plt.plot(unstab_up['I_app'],unstab_up['V'],'b--')
plt.plot(unstab_down['I_app'],unstab_down['V'],'b--')
plt.plot(vmax['I_app'],vmax['V'],'b-',label='$V_{max}$')
plt.plot(vmin['I_app'],vmin['V'],'b-',label='$V_{min}$')
plt.vlines(x=stable[-1][0],ymin=min_v[0][1],ymax=max_v[0][1],linestyles=':',color='b')
plt.text(90,-65,'$V_{min}$'); plt.text(90,20,'$V_{max}$')
plt.ylim(-90.,40.)
plt.xlabel('$I_{app}$',fontsize=18); plt.ylabel('V',fontsize=18)
plt.grid(linestyle=':')
plt.show()

# %%