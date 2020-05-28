# This script calculates the delta-v based on these transfers:
# - circularisation
# - 80km corrections
# - phase shift towards constellation
# - realignment
# - orbit maintenances
# - End Of Life manoeuvre

import numpy as np

# Orbital parameters
m_earth = 5.9722 * 10 ** 24                             # kg
G = 6.67 * 10 ** -11                                    #
mu_earth = m_earth * G / (10 ** 9)                      #
r_earth = 6378.136                                      # km
r_GEO = 42164.136                                       # km
T_GEO = 2 * np.pi * np.sqrt(r_GEO ** 3 / mu_earth)      # s
omega_GEO = np.sqrt(mu_earth / (r_GEO ** 3))            # rad/s
peri_GTO = 250      # altitude                          # km
apo_GTO = 35786     # altitude                          # km
a_GTO = (peri_GTO + apo_GTO) / 2 + r_earth              # km
T_GTO = 2 * np.pi * np.sqrt(a_GTO ** 3 / mu_earth)      # s
d = np.sqrt(r_GEO ** 2 * 3)                             # km


def Hohmann(r_dep, r_tar, mu_earth, r_earth):
    # This function calculates the total delta-v needed for a Hohmann transfer,
    # knowing the initial orbit and the target orbit (r_dep and r_tar)
    a_T = 0.5*(r_dep+r_tar)
    Vc_dep = np.sqrt(mu_earth/(r_dep + r_earth))
    Vc_tar = np.sqrt(mu_earth/(r_tar + r_earth))
    dV_1 = np.sqrt(mu_earth*(2/(r_dep + r_earth) - 1/(a_T + r_earth)))-Vc_dep
    dV_2 = np.sqrt(mu_earth*(2/(r_tar + r_earth) - 1/(a_T + r_earth)))-Vc_tar

    return abs(dV_1) + abs(dV_2)


def Phase_Shift(shift, t_transfer, r_GEO, T_GEO, mu_earth):
    dT = t_transfer * T_GEO
    #dpos = shift / (2 * np.pi * r_GEO) * 2 * np.pi
    da = shift * r_GEO / (3 * np.pi) * T_GEO / dT
    dV = 2 * (np.sqrt(mu_earth * (2 / r_GEO - 1 / (r_GEO + da))) - np.sqrt(mu_earth / r_GEO))
    return dT, da, dV

# circularisation


# 80 km corrections
r_dep, r_tar = 35706, 35786
eighty_km = Hohmann(r_dep, r_tar, mu_earth, r_earth)
print(eighty_km)
circularisation_orbit = Hohmann(r_dep=1, r_tar=2, mu_earth=398600, r_earth=6378)


# phase shift towards 120 degrees apart
shift_120 = 2/3 * np.pi - (3 * omega_GEO * T_GTO - 2 * np.pi)
t_transfer_120 = 7
dT_120, da_120, dV_120 = Phase_Shift(shift_120, t_transfer_120, r_GEO, T_GEO, mu_earth)

# realignment

# orbit maintenance

# end-of-life manoeuvres
