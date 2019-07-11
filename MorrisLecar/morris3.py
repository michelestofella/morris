# =============================================================================
# 
# SCRIPT FOR BIFURCATION ANALYSIS OF THE MORRIS LECAR MODEL
# Runge Kutta algorithm is applied for different I_app values
# Frequency of the output signal is calculated
# Frequency plot is generated
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

""" Set I_app values and v0_values to draw bifurcation diagram """
I_app_values = np.linspace(0,100,51)    

""" Integrate the model at different values of I_app 
    to find the frequency of the generated action potential """
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
print("Done: 100.0 %")

# %%

""" Plot of frequencies vs I_app """
plt.plot(I_app_values,freq)
plt.xlim(0,100); plt.ylim(0,160)
plt.xlabel('$I_{app}$', fontsize=18)
plt.ylabel('Frequency [Hz]', fontsize=18)
plt.grid(linestyle=':')

# %%