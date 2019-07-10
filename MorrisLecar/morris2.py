""" Import Libraries """
import numpy as np
import matplotlib.pyplot as plt  
import pandas as pd
from scipy.signal import find_peaks

""" Import algorithms """
import sys
sys.path.insert(0, '../newton')
from newton2 import newton2
sys.path.insert(0, '../rungekutta')
from rk4_system import RK4_system

""" Fixed Model Parameters """
g_ca = 20.0;     g_k = 20.;   g_leak = 2.
E_ca = 50.0;    E_k = -100.; E_leak = -70.
phi_w = 0.15
c = 2.

v_ca  = -12.0;   v_k = -10.
theta_ca = 18.0; theta_k = 13.

# %%

""" Model for the Runge Kutta Integration """
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

# %%

""" Model for Bifurcation Analysis """
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
I_app_values = np.linspace(0,100,201)    
v0_values = np.linspace(-80,40,61)
w0 = 0.0

""" Find the zeros of the function """
v_zeros = []; w_zeros = []
for I_app in I_app_values:
    zero_v = []; zero_w = []
    zeros_v = []; zeros_w = []
    for v0 in v0_values:
        p0 = [v0,w0]
        point = newton2(f,Jf,p0)
        if point is not False:
            zero_v.append(round(point[0],5))
            zero_w.append(round(point[1],5))
        for element in zero_v:
            if element not in zeros_v:
                zeros_v.append(element)
        for element in zero_w:
            if element not in zeros_w:
                zeros_w.append(element)
    v_zeros.append(zeros_v); w_zeros.append(zeros_w)
print('Several Errors may be shown: do not worry about them,')
print('they are errors encountered by the Newton algorithm that returns no values.')

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
    to find maximum and minimum values of the limit cycle
    and to find the frequency of the generated action potential """
g = [g1,g2]; dt = 0.01; t0 = 0.0; Nstep = 1000
y0 = [-25.0,0.0]

max_v = []; min_v = []; freq = []
for j in range(0,len(I_app_values)):
    I_app = I_app_values[j]
    print("Calculating:",round(j*100/len(I_app_values),1),"%")
    time, sol = RK4_system(g, dt, y0, t0, Nstep)
    peaks, _ = find_peaks(sol[0], height=0)
    peaks = peaks*dt/1000
    period = []
    for i in range(1,len(peaks)):
        period.append(peaks[i]-peaks[i-1])
    freq.append(1/np.mean(period))
    if I_app > stable[-1][0]:
        max_v.append([I_app,max(sol[0])])
        min_v.append([I_app,min(sol[0])])
print("Done: 100.0 %")
                
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

# %%

""" Integrate the model for two different situations """
g = [g1,g2]; dt = 0.01; t0 = 0.0; Nstep = 5000

y0_1 = [-25.0,0.0]
y0_2 = [-10.0,-0.0]
I_app = 14.
time, sol = RK4_system(g, dt, y0_1, t0, Nstep)
time2,sol2 = RK4_system(g, dt, y0_2, t0, Nstep)

# %%

""" Plot of signal and phase space """
delta = 0.025
x = np.arange(-100.0, 50.0, delta)
y = np.arange(-0.1, .4, delta)
X, Y = np.meshgrid(x, y)
Z1 = (g_ca*(0.5*(1+np.tanh((X-v_ca)/theta_ca)))*(E_ca-X) + g_k*Y*(E_k-X) +\
     g_leak*(E_leak-X) + I_app)/c
Z2 = phi_w*((0.5*(1+np.tanh((X-v_k)/theta_k))-Y)/(1/(np.cosh((X-v_k)/(2*theta_k)))))

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(15,6))
ax1.plot(sol[0],sol[1], color='blue')
ax1.plot(sol2[0],sol2[1], color='red')
ax1.contour(X, Y, Z1, 0, colors='black', linestyles='--', linewidths=1)
ax1.contour(X, Y, Z2, 0, colors='black', linestyles='--', linewidths=2)
ax1.set_xlabel('V',fontsize=18); ax1.set_ylabel('w',fontsize=18)
ax1.set_ylim([-0.1,0.4])
ax2.plot(time, sol[0], 'b')
ax2.plot(time2,sol2[0],'r')
ax2.set_xlabel('Time', fontsize=18); ax2.set_ylabel('Voltage',fontsize=18)
ax1.grid(linestyle=':'); ax2.grid(linestyle=':')

# %%

""" Plot of frequencies vs I_app """
plt.plot(I_app_values,freq)
plt.xlim(0,100); plt.ylim(0,160)
plt.xlabel('$I_{app}$', fontsize=18)
plt.ylabel('Frequency [Hz]', fontsize=18)
plt.grid(linestyle=':')

# %%