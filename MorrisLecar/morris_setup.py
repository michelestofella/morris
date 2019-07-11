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
g_ca = 20.0;     g_k = 20.;   g_leak = 2.
E_ca = 50.0;    E_k = -100.; E_leak = -70.
phi_w = 0.15
c = 2.
v_k = -10.; theta_ca = 18.0; theta_k = 13.

''' Analysis is performed in the (I_app,v_ca) plane '''
v_ca  = 0.0  

# %%