# =============================================================================
# 
# SCRIPT FOR NUMERICAL INTEGRATION OF MORRIS LECAR MODEL
# Runge Kutta algorithm is applied 
# Plots of the signal and the phase space are generated
#
# =============================================================================

from morris_setup import *

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

# %%

""" Integrate the model for two different situations """
g = [g1,g2]; dt = 0.01; t0 = 0.0; Nstep = 5000

y0_1 = [-25.0,0.0]
y0_2 = [-10.0,-0.0]
I_app = 14.0

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