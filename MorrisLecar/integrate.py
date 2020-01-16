# =============================================================================
# 
# SCRIPT FOR NUMERICAL INTEGRATION OF MORRIS LECAR MODEL
# Runge Kutta algorithm is applied 
# Plots of the signal and the phase space are generated
#
# Several parameters must be parsed in order to run the algorithm:
#   --v_ca  the parameter of the morris lecar model that discriminates
#           different classes of neurons.
#   --I_app the external current applied to the model 
#   --dt    integration time step
#   --Nstep number of integration steps to be performed
#   --v0    initial condition on the voltage
#   --w0    initial condition of the fraction of opened channels
#   --out   name of the generated figure
#
# =============================================================================

import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--v_ca")
parser.add_argument("--I_app")
parser.add_argument("--dt")
parser.add_argument("--Nstep")
parser.add_argument("--v0")
parser.add_argument("--w0")
parser.add_argument("--out")

config = {}
opts = parser.parse_args()

if opts.v_ca:
    v_ca = float(opts.v_ca)
if opts.I_app:
    I_app = float(opts.I_app)
if opts.dt:
    dt = float(opts.dt)
else:
    dt = 0.01
if opts.Nstep:
    Nstep = int(opts.Nstep)
else:
    Nstep = 10000
if opts.v0:
    v0 = float(opts.v0)
else:
    v0 = 0
if opts.w0:
    w0 = float(opts.w0)
else:
    w0 = 0
if opts.out:
    save = True
    out = opts.out
else:
    save = False

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

# %%

""" Integrate the model """
g = [g1,g2]
t0 = 0.0
y0 = [v0,w0]

time, sol = RK4_system(g, dt, y0, t0, Nstep)

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
ax1.contour(X, Y, Z1, 0, colors='black', linestyles='--', linewidths=1)
ax1.contour(X, Y, Z2, 0, colors='black', linestyles='--', linewidths=2)
ax1.set_xlabel('V',fontsize=18); ax1.set_ylabel('w',fontsize=18)
ax1.set_ylim([-0.1,0.4])
ax2.plot(time, sol[0], 'b')
ax2.set_xlabel('Time', fontsize=18); ax2.set_ylabel('Voltage',fontsize=18)
ax1.grid(linestyle=':'); ax2.grid(linestyle=':')
if save == True:
    plt.savefig(out+'.png')
plt.show()

# %%