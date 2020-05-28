# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np

# Chemical propulsion (impulsive manoeuvres)
# - circularisation
# - phase shift towards constellation
# - realignment
# - other manoeuvres that are not possible with electric

# Electric propulsion (non-impulsive manoeuvres)
# - 80km corrections?
# - realignment?
# - orbit maintenances?
# - End Of Life manoeuvre?

# time for transfer
def Transfer_Time(m, I_sp, T, r, r0)
t = m * g * I_sp / T * (1 - np.exp(1 / (I_sp * g)))