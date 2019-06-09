# =============================================================================
#
# MORRIS LECAR MODEL
#
# Algorithm tests:      !pytest rk4_system.py
# 
# Source for the model: LIU, 'Bifurcation Analysis of a Morris-Lecar Neuron 
#                       Model', Biological Cybernetics (2014)
#
# =============================================================================

import numpy as np
import matplotlib.pyplot as plt    

""" Runge-Kutta algorithm """
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

# %%

""" Morris-Lecar Model """
def f1(t,y1,y2):
    part1 = g_ca*m_inf(y1)*(E_ca-y1)
    part2 = g_k*y2*(E_k-y1)
    part3 = g_leak*(E_leak-y1)
    return (part1+part2+part3+I_app)/c
def f2(t,y1,y2):
    return phi_w*((w_inf(y1)-y2)/tau_w(y1))

def m_inf(v):
    return 0.5*(1+np.tanh((v-v_ca)/theta_ca))
def w_inf(v):
    return 0.5*(1+np.tanh((v-v_k)/theta_k))
def tau_w(v):
    return 1/(np.cosh((v-v_k)/(2*theta_k)))
    
""" Model Parameter """
g_ca = 20.0;     g_k = 20.;   g_leak = 2.
E_ca = 50.0;    E_k = -100.; E_leak = -70.
phi_w = 0.15
c = 2.

v_ca  = 0.;   v_k = -10.
theta_ca = 18.0; theta_k = 13.
I_app = 80.

""" Integration Parameters """
f = [f1,f2]; dt = 0.01; t0 = 0.0; Nstep = 5000
y0 = [-25.0,0.0]

""" Application of the algorithm """
time, y = RK4_system(f, dt, y0, t0, Nstep)

# %%

""" Data visualization """
f, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15,5))
""" Membrane voltage vs Time """
ax1.plot(time, y[0])
""" Potassium gating variable vs Time """
ax2.plot(time, y[1])
""" Potassium gating variable w vs Membrane Voltage V """
ax3.plot(y[0], y[1])

# %%