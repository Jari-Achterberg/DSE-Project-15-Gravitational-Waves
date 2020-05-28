# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np
from Astrodynamics import mu_earth, dV_EOL, dV_eighty_km, dV_circularisation, dV_maintenance, dV_120, dV_realignment, total_delta_v

# input parameters
m_i = 20    # kg
g0 = 9.81


def Transfer_Time(m, I_sp, T, r, r0):
    # time for transfer
    g = g0
    t = m * g * I_sp / T * (1 - np.exp(1 / (I_sp * g) * (np.sqrt(mu_earth/r) - np.sqrt(mu_earth/r0))))
    return t


def rocket_equation(delta_v, isp, m_i):

    w = g0*isp
    Mp = (1 - np.exp(-delta_v/w))*m_i
    return Mp


# Chemical propulsion (impulsive manoeuvres)
# Only circularisation

delta_v = dV_circularisation * 1000
delta_v_2 = total_delta_v * 1000
# all options considered from left to right: (6 options)
# Hydros-C, Lunar Flashlight MiPS, Argomoon Hibrid MiPS, Green Hybdrid, GR-1, 1 N GPHP, EPSS 1CKK
isps = [310, 190, 190, 215, 231, 204, 252]
Mp_list, Mp_list_2 = [], []
for isp in isps:
    Mp = rocket_equation(delta_v, isp, m_i)
    Mp_2 = rocket_equation(delta_v_2, isp, m_i)
    Mp_list.append(Mp)
    Mp_list_2.append(Mp_2)
print(Mp_list)
print(Mp_list_2)
# All manoeuvres

# Electric propulsion (non-impulsive manoeuvres)
# - 80km corrections?
# - phase shift towards constellation?
# - realignment?
# - orbit maintenances
# - End Of Life manoeuvre?

m = 20



