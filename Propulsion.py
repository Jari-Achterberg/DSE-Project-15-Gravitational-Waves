# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np
from matplotlib import pyplot as plt
from Astrodynamics import mu_earth, dV_EOL, dV_eighty_km, dV_circularisation, dV_maintenance, dV_120, dV_realignment, total_delta_v

# input parameters
m_i = 20    # kg
g0 = 9.81


def Transfer_Time(m, I_sp, T, r, r0):
    # time for transfer
    g = g0 / 1000
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
# Hydros-C, Lunar Flashlight MiPS, Argomoon Hibrid MiPS, Green Hybdrid, GR-1, 1 N GPHP, EPSS 1CKK # chemical
isps = [310, 190, 190, 215, 231, 204, 252, 3000]
Mp_list, Mp_list_2 = [], []
for isp in isps:
    Mp = rocket_equation(delta_v, isp, m_i)
    Mp_2 = rocket_equation(delta_v_2, isp, m_i)
    Mp_list.append(Mp)
    Mp_list_2.append(Mp_2)

print(Mp_list)
print(Mp_list_2)
# All manoeuvres
T = [1.2, 0.4, 0.1, 8, 0.26, 0.25, 1, 50/(10**6)]
# Electric propulsion (non-impulsive manoeuvres)
# - 80km corrections?
# - phase shift towards constellation?
# - realignment?
# - orbit maintenances
# - End Of Life manoeuvre?

m = m_i
# PPTCUP, mu-CAT, PUC, Resistojet System, electrospray, MPACS,
isps_2 = [655, 3000, 70, 99, 800, 827]
T_2 = [40/(10**6), 50/(10**6), 0.45, 100/1000, 0.7/1000, 80/(10**6)]
r = 42164
r0 = 250 + 6378

t_list = []

for i in range(len(isps_2)):
    t = Transfer_Time(m, isps_2[i], T_2[i] / 1000, r, r0)
    t /= 60
    t /= 60
    t /= 24
    t_list.append(t)


print("transfer times: ", t_list)

Isp = 100
T_test = 0
t_test_list = []
T_test_list = []
for i in range(100):
    T_test += 1/10000
    t_test = Transfer_Time(m, Isp, T_test / 1000, r, r0)
    t_test /= 3600
    t_test /= 24
    T_test_list.append(T_test)
    t_test_list.append(t_test)

print(T_test_list[-1], t_test_list[-1])
plt.plot(T_test_list, t_test_list)
plt.show()
