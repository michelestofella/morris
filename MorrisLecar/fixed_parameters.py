""" Import Libraries """
import numpy as np
import matplotlib.pyplot as plt  
import pandas as pd
from scipy.signal import find_peaks

""" Import algorithms """
import sys
sys.path.insert(0, '../rungekutta')
from rk4_system import RK4_system
sys.path.insert(1, '../newton')
from newton2 import newton2

""" Fixed Model Parameters """
# calcium channels
g_ca = 20.0         # maximal conductance
E_ca = 50.0         # Nernst potential
theta_ca = 18.0     # steady state parameter
# potassium channels
g_k = 20.           # maximal conductance
E_k = -100.         # Nernst potential
v_k = -10.          # steady state parameter
theta_k = 13.       # steady state parameter
# leak current 
g_leak = 2.         # maximal conductance
E_leak = -70.       # Nernst potential
# other parameters
c = 2.              # membrane conductance
phi_w = 0.15        # constant term 

# %%