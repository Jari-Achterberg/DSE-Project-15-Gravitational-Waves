# This script calculates the delta-v based on these transfers:
# - circularisation
# - 80km corrections
# - phase shift towards constellation
# - realignment
# - orbit maintenances
# - End Of Life manoeuvre

import numpy as np

# Orbital parameters
m_earth = 5.9722 * 10 ** 24
G = 6.67 * 10 ** -11
mu_earth = m_earth * G / (10 ** 9)
r_earth = 6378.136
r_GEO = 42164.136
T_GEO = 2 * np.pi * np.sqrt(r_GEO ** 3 / mu_earth)
omega_GEO = np.sqrt(mu_earth / (r_GEO ** 3))
peri_GTO = 250 #altitude
apo_GTO = 35786 #altitude
a_GTO = (peri_GTO + apo_GTO) / 2 + r_earth
T_GTO = 2 * np.pi * np.sqrt(a_GTO ** 3 / mu_earth)
d = np.sqrt(r_GEO ** 2 * 3)




def Hohmann():
    t = 0
    return t


def Phase_Shift(shift, t_transfer, r_GEO, T_GEO, mu_earth):
    dT = t_transfer * T_GEO
    dpos = shift / (2 * np.pi * r_GEO) * 2 * np.pi
    da = dpos * r_GEO / (3 * np.pi) * T_GEO / dT
    dV = 2 * (np.sqrt(mu_earth * (2 / r_GEO - 1 / (r_GEO + da))) - np.sqrt(mu_earth / r_GEO))
    return dV


circularisation_orbit = Hohmann()

# phase shift towards 120 degrees apart
shift_120 = 2/3 * np.pi - (3 * )
t_transfer_120 = 7
dV_realignment
