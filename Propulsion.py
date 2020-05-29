# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np
from matplotlib import pyplot as plt
from Astrodynamics import mu_earth, dV_EOL, dV_eighty_km, dV_circularisation, dV_maintenance, dV_120, dV_realignment, total_delta_v
from scipy import integrate

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


def Time_Transfer():
    # assume e0 for now
    e0 = 0.9897524756
    omega = 178 / 180*np.pi
    af = 42164.136  # km
    i0 = 6.02 / 180*np.pi
    TT = 0.3/1000
    m = 20
    a0 = 24396  # km

    v = 1.64
    r = af
    F_in = m*v**2/r
    F_min = 0.00029774
    k = F_min/F_in
    print(F_min)
    print(F_in)
    F = TT/m

    getal = np.sqrt((1 - k) / 2)
    print(getal)
    K = integrate.quad(lambda t: 1/(np.sqrt(1-getal**2*t**2)*np.sqrt(1-t**2)), 0, 1)
    E = integrate.quad(lambda t: np.sqrt(1-getal**2*t**2)/np.sqrt(1-t**2), 0, 1)

    ksi = 4/3 * np.sqrt(2/(1-k)) * ((1+k)*K[0]-2*k*E[0])

    u1 = np.pi*(1+k)/(2*ksi)
    u2 = np.pi*(7+5*k)/(12*ksi)

    a = lambda x: af * np.exp(2/(1+9*u2**2)*((3*u2-2*u1)*np.arcsin(x)-(1+6*u1*u2)*np.log(3*u2*x+np.sqrt(1-x**2))))  # a0/af
    print("a0: ", a(e0))
    H = integrate.quad(lambda x: 1 / (np.sqrt(a(x)) * ((1-x**2)+3*u2*x*np.sqrt(1-x**2))), 0, e0)                     # H
    fc = 1/(6*u2*(1+9*u2))*((18*u2**2+4)*np.log(np.abs(np.sqrt(1-e0**2)+3*u2*e0))-2*(1+9*u2**2)*np.log(1-e0**2)-6*u2*np.arcsin(e0))
    fs = 1/(np.sqrt(1+9*u2**2))*(2*np.log(np.abs(np.sqrt(1+9*u2**2)+3*u2))-np.log(np.abs(((np.sqrt(1+9*u2**2)+3*u2)*(1+np.sqrt(1-e0**2))-3*u2*e0)/((np.sqrt(1+9*u2**2)+3*u2)*(1+np.sqrt(1-e0**2))+3*u2*e0))))

    F_in = F * (1 + (i0**2*ksi**2)/(16*(fc*np.cos(omega)**2+fs*np.sin(omega)**2)**2)) ** (-1/2)
    t = 2 * np.pi * np.sqrt(mu_earth) * H[0] / (F_in * ksi)
    print(H)
    # dv = 2*np.pi*H[0]/ksi*np.sqrt(mu_earth+(mu_earth*i0**2*ksi**2)/(16*(fc*np.cos(omega)**2+fs*np.sin(omega)**2)**2))
    # t = dv/F
    return t


tt = Time_Transfer()
print(tt)
'''
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
'''
