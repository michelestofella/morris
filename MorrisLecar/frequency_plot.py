# =============================================================================
# 
# SCRIPT FOR BIFURCATION ANALYSIS OF THE MORRIS LECAR MODEL
# Runge Kutta algorithm is applied for different I_app values
# Frequency of the output signal is calculated
# Frequency plot is generated
# 
# Several parameters must be parsed in order to run the algorithm:
#   --v_ca  the parameter of the morris lecar model that discriminates
#           different classes of neurons.
#   --dt    integration time step
#   --Nstep number of integration steps to be performed
#   --v0    initial condition on the voltage
#   --w0    initial condition of the fraction of opened channels
#   --Imin  minimum value of applied current analysed
#   --Imax  maximum value of applied current analysed 
#   --out   name of the generated figure
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
parser.add_argument("--out")

config = {}
opts = parser.parse_args()

if opts.v_ca:
    v_ca = float(opts.v_ca)
if opts.dt:
    dt = float(opts.dt)
else:
    dt = 0.01
if opts.Nstep:
    Nstep = int(opts.Nstep)
else:
    Nstep = 5000
if opts.v0:
    v0 = float(opts.v0)
else:
    v0 = 0
if opts.w0:
    w0 = float(opts.w0)
else:
    w0 = 0

if opts.Imin:
    Imin = float(opts.Imin)
if opts.Imax:
    Imax = float(opts.Imax)
    
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

""" Set I_app values and v0_values to draw bifurcation diagram """
I_app_values = np.linspace(Imin,Imax,2*(Imax-Imin)+1)    

""" Integrate the model at different values of I_app 
    to find the frequency of the generated action potential """
g = [g1,g2]
t0 = 0.0
y0 = [v0,w0]

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
print("Done: 100.0 %")

# %%

""" Plot of frequencies vs I_app """
plt.figure(figsize=(15,10))
plt.plot(I_app_values,freq)
plt.xlim(0,100); plt.ylim(0,160)
plt.xlabel('$I_{app}$', fontsize=18)
plt.ylabel('Frequency [Hz]', fontsize=18)
plt.grid(linestyle=':')
if save == True:
    plt.savefig(out+'.png')
plt.show()

# %%