# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np
from Astrodynamics import mu_earth, dV_EOL, dV_eighty_km, dV_circularisation, dV_maintenance, dV_120, dV_realignment

# Chemical propulsion (impulsive manoeuvres)
# Only circularisation


# All manoeuvres



# Electric propulsion (non-impulsive manoeuvres)
# - 80km corrections?
# - phase shift towards constellation?
# - realignment?
# - orbit maintenances
# - End Of Life manoeuvre?

m = 20


def Transfer_Time(m, I_sp, T, r, r0):
    # time for transfer
    g = 9.81
    t = m * g * I_sp / T * (1 - np.exp(1 / (I_sp * g) * (np.sqrt(mu_earth/r) - np.sqrt(mu_earth/r0))))
    return t


def rocket_equation(delta_v, isp, m_i):
    g0 = 9.81
    w = g0*isp
    Mp = (1 - np.exp(-delta_v/w))*m_i
    return Mp

