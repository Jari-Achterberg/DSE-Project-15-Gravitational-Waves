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
V_apo_GTO = 1.64                                        # km/s
i_GTO = 6.02                                            # degrees

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
Vc_GEO = np.sqrt(mu_earth/r_GEO)
dV_circularisation = np.sqrt(V_apo_GTO ** 2 + Vc_GEO ** 2 - 2 * V_apo_GTO * Vc_GEO * np.cos(i_GTO * np.pi/180))

# 80 km corrections
r_dep, r_tar = 35706, 35786
dV_eighty_km = Hohmann(r_dep, r_tar, mu_earth, r_earth)


# phase shift towards 120 degrees apart
# shift_120 = 2/3 * np.pi - (3 * omega_GEO * T_GTO - 2 * np.pi)
shift_120 = -2/3 * np.pi
t_transfer_120 = 7
dT_120, da_120, dV_120 = Phase_Shift(shift_120, t_transfer_120, r_GEO, T_GEO, mu_earth)

# realignment
# only hohmann is considered, as it results in the highest delta v
max_drift = np.sin(5 / 60 * np.pi / 180)*d
r1_rea = apo_GTO - max_drift
r2_rea = apo_GTO
dV_rea = Hohmann(r1_rea, r2_rea, mu_earth, r_earth)

number_of_realignments = 49  # based on three years
dV_realignment = dV_rea*number_of_realignments

# orbit maintenance
dV_orbit = 2 * 0.0075
dV_attitude = 2 * 0.006
dV_momentum = 2 * 0.006
dV_maintenance = dV_orbit + dV_attitude + dV_momentum

# end-of-life manoeuvres
r1_EOL = apo_GTO            # == GEO
r2_EOL = r1_EOL + 300       # assume EOL orbit is 300 km higher
dV_EOL = Hohmann(r1_EOL, r2_EOL, mu_earth, r_earth)
# TEST
#### total, correction of 10% is assumed
total_delta_v = (dV_EOL + dV_eighty_km + dV_realignment + dV_maintenance + abs(dV_120) + dV_circularisation)*1.1
dV_in_orbit = total_delta_v - 1.1*dV_circularisation
