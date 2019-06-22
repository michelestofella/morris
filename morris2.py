""" Import Libraries """
import numpy as np
import matplotlib.pyplot as plt  
from numpy.linalg import inv
import pandas as pd
  
""" Fixed Model Parameters """
g_ca = 20.0;     g_k = 20.;   g_leak = 2.
E_ca = 50.0;    E_k = -100.; E_leak = -70.
phi_w = 0.15
c = 2.

v_ca  = 0.0;   v_k = -10.
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
    
""" Runge Kutta Algorithm """
def RK4_system(f, dt, y0, t0, Nstep):
    Neq = len(f)
    t = np.zeros(Nstep+1); t[0] = t0
    y = np.zeros([Neq, Nstep+1]) 
    dy1 = np.zeros([Neq, Nstep+1]); dy2 = np.zeros([Neq, Nstep+1]); 
    dy3 = np.zeros([Neq, Nstep+1]); dy4 = np.zeros([Neq, Nstep+1])
    for j in range(0, Neq):
        y[j][0] = y0[j]
    
    for i in range(0,Nstep):
        for j in range(0,Neq):
            dy1[j] = dt*f[j](t,*y)
            dy2[j] = dt*f[j](t+0.5*dt,*y+0.5*dy1)
            dy3[j] = dt*f[j](t+0.5*dt,*y+0.5*dy2)
            dy4[j] = dt*f[j](t+dt,*y+dy3)
        
            t[i+1] = t[i] + dt
            y[j][i+1] = y[j][i]+(dy1[j][i]+2*dy2[j][i]+2*dy3[j][i]+dy4[j][i])/6
    return t, y

""" Newton Algorithm """    
def newton2(f,Jf,p0,eps=1e-8,max_iter=20):
    for k in range(0,max_iter):
        f_k = f(p0[0],p0[1])
        Jf_k = Jf(p0[0],p0[1])
        invJf_k = inv(Jf_k)
        pk = p0 - np.dot(invJf_k,f_k)
    
        dist = np.sqrt(np.dot(pk-p0,pk-p0))
        if dist < eps:
            return pk
        p0 = pk
    return False

# %%

I_app_values = np.linspace(0,100,100)    
v0_values = np.linspace(-80,40,100)
w0 = 0.

zero_v = []; zero_w = [];
zeros_v = []; zeros_w = []
for I_app in I_app_values:
    for v0 in v0_values:
        p0 = [v0,w0]
        point = newton2(f,Jf,p0)
        if point is not False:
            zero_v.append(round(point[0],6))
            zero_w.append(round(point[1],6))
        for element in zero_v:
            if element not in zeros_v:
                zeros_v.append(element)
        for element in zero_w:
            if element not in zeros_w:
                zeros_w.append(element)
print('Several Errors may be shown: do not worry about them, they are errors encountered by the Newton algorithm that return no values.')

stable = []; unstable = []
for i in range(0,len(I_app_values)):
    I_app = I_app_values[i]
    eigen = np.linalg.eigvals(Jf(zeros_v[i],zeros_w[i]))
    for j in range(0,len(eigen)):
        if eigen[j]<0:
            stable.append([I_app_values[i],zeros_v[i],zeros_w[i],'stable'])
        else:
            unstable.append([I_app_values[i],zeros_v[i],zeros_w[i],'unstable'])
            
# %%

""" Integration Parameters """
g = [g1,g2]; dt = 0.01; t0 = 0.0; Nstep = 5000
y0 = [-25.0,0.]

max_v = []; min_v = []
for I_app in I_app_values:
    if I_app > unstable[0][0]:
        time, sol = RK4_system(g, dt, y0, t0, Nstep)
        max_v.append([I_app,max(sol[0])])
        min_v.append([I_app,min(sol[0])])
        
# %%            
 
col_st = ['I_app','V','w','Stability']
stab = pd.DataFrame(stable,columns=col_st)
unstab = pd.DataFrame(unstable, columns=col_st)

col=['I_app','V']
vmax = pd.DataFrame(max_v,columns=col)
vmin = pd.DataFrame(min_v,columns=col)

plt.plot(stab['I_app'],stab['V'])
plt.plot(unstab['I_app'],unstab['V'])
plt.plot(vmax['I_app'],vmax['V'],'b-',label='$V_{max}$')
plt.plot(vmin['I_app'],vmin['V'],'b-',label='$V_{min}$')
plt.xlabel('$I_{app}$',fontsize=18); plt.ylabel('V',fontsize=18)
plt.vlines(x=unstable[0][0],ymin=min_v[0][1],ymax=max_v[0][1],linestyles=':',color='b')
plt.grid(linestyle=':')

# %%