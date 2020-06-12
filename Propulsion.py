# Determine what type of propulsion system is suitable for which manoeuvre

import numpy as np
from matplotlib import pyplot as plt
from Astrodynamics import mu_earth, total_delta_v, dV_in_orbit, dV_circularisation
from scipy import integrate


def Transfer_Time(m, I_sp, T, r, r0):
    # this function calculates the time it takes for transfer from one circular orbit to another circular orbit
    # still not sure if it is even correct (to be verified)

    g = g0 / 1000
    t = m * g * I_sp / T * (1 - np.exp(1 / (I_sp * g) * (np.sqrt(mu_earth/r) - np.sqrt(mu_earth/r0))))
    return t


def rocket_equation(isp, m_i, v, rho):
    mp = (rho * v) / 1000

    w = isp * g0
    m_f = m_i - mp
    print(m_i, m_f)
    delta_v2 = w * np.log(m_i/m_f)
    return delta_v2


def rocket_equation_reversed(delta_v, isp, m_i):
    # calculation of propellant mass using the rocket equation with inputs:
    # - velocity difference (delta v), specific impulse (isp) and initial mass (m_i)
    w = g0*isp
    Mp = (1 - np.exp(-delta_v/w))*m_i
    return Mp


def Time_Transfer():
    # NOT USED ANYMORE
    # This function calculates the time it takes for a transfer from an arbitrary orbit to GEO
    # inputs are given within this function, should be adjusted in the future
    # inputs are: eccentricity of initial orbit (e0), semi major axis initial and final orbit (a0 and af),
    # inclination of initial orbit (i0), thrust force of electrical system (TT), initial mass of s/c (m)
    # k is found through iteration
    # -----------------------------------------------------------------------------------------------------
    # latest changes: (01-06)
    # - changed e0
    # - changed thrust
    # - changed k based on e0 difference
    # - changed np.log to np.log10 because log is used in function instead of ln
    # - changed k accordingly

    # test uses values from paper
    test = False
    if test:
        e0 = 0.44
        omega = 90 / 180 * np.pi
        i0 = 90 / 180 * np.pi
        TT = 1*10**(-5)
        m = 50

    else:
        # e0 = 0.9897524756
        e0 = 0.73
        omega = 178 / 180 * np.pi
        i0 = 6.02 / 180 * np.pi
        TT = 100/1000/1000
        m = 20

    # a0 = 24396 km
    af = 42164.136  # km

    # F_in = TT/m
    # k = F_min/F_in
    # k = 0.233379496         # for e0 = 0.989
    # k = 0.02438             # for e0 = 0.73
    k = 0.69682                # for e0 = 0.73 AND np.log change
    # k = 0.987
    F = TT/m

    # functions K and E
    nr = np.sqrt((1 - k) / 2)
    K = integrate.quad(lambda t1: 1/(np.sqrt(1-nr**2*t1**2)*np.sqrt(1-t1**2)), 0, 1)
    E = integrate.quad(lambda t1: np.sqrt(1-nr**2*t1**2)/np.sqrt(1-t1**2), 0, 1)

    # ksi, u1, u2,
    ksi = 4/3 * np.sqrt(2/(1-k)) * ((1+k)*K[0]-2*k*E[0])
    u1 = np.pi*(1+k)/(2*ksi)
    u2 = np.pi*(7+5*k)/(12*ksi)

    # a0/af
    a = lambda x: af * np.exp(2/(1+9*u2**2)*((3*u2-2*u1)*np.arcsin(x)-(1+6*u1*u2)*np.log10(3*u2*x+np.sqrt(1-x**2))))
    print("a0: ", a(e0))

    # H / fc / fs
    H = integrate.quad(lambda x: 1 / (np.sqrt(a(x)) * ((1-x**2)+3*u2*x*np.sqrt(1-x**2))), 0, e0)
    fc = 1/(6*u2*(1+9*u2))*((18*u2**2+4)*np.log(np.abs(np.sqrt(1-e0**2)+3*u2*e0))-2*(1+9*u2**2)*np.log(1-e0**2)-6*u2*np.arcsin(e0))
    fs = 1/(np.sqrt(1+9*u2**2))*(2*np.log(np.abs(np.sqrt(1+9*u2**2)+3*u2))-np.log(np.abs(((np.sqrt(1+9*u2**2)+3*u2)*(1+np.sqrt(1-e0**2))-3*u2*e0)/((np.sqrt(1+9*u2**2)+3*u2)*(1+np.sqrt(1-e0**2))+3*u2*e0))))
    print(fc)
    print(fs)
    # F_in = F * (1 + (i0**2*ksi**2)/(16*(fc*np.cos(omega)**2+fs*np.sin(omega)**2)**2)) ** (-1/2)
    # t = 2 * np.pi * np.sqrt(mu_earth) * H[0] / (F_in * ksi)

    # dv and t
    dv = 2*np.pi*H[0]/ksi*np.sqrt(mu_earth)*np.sqrt(1+(i0**2*ksi**2)/(16*(fc*(np.cos(omega))**2+fs*(np.sin(omega))**2)**2))
    print("part: ", (i0**2*ksi**2)/(16*(fc*(np.cos(omega))**2+fs*(np.sin(omega))**2)**2))
    print("dv: ", dv)
    t = dv/F
    return t


# input parameters
m_i = 16.5    # kg
g0 = 9.81


# Chemical propulsion (impulsive manoeuvres)
# delta_v == only
# check for usage of including transfer (delta-v-2) or excluding transfer (delta-v-1)
delta_v = dV_in_orbit * 1000
delta_v_2 = total_delta_v * 1000
# all options considered from left to right: (6 options)
# Hydros-C, Lunar Flashlight MiPS, Argomoon Hibrid MiPS, Green Hybdrid, GR-1, 1 N GPHP, EPSS C1K    # chemical
# specific impulses, thrust levels and densities per option:
isps = [310, 190, 190, 215, 231, 204, 213]
rho = [1.0, 1.1, 1.1, 1.1, 1.1, 1.2648, 1.2648]
T = [1.2, 0.4, 0.1, 8, 0.26, 0.25, 0.22, 50/(10**6)]

# calculate propellant mass of each propellant type for the two delta-v options
# calculate burn time with the found propellant mass
Mp_list, Mp_list_2 = [], []
for i, isp in enumerate(isps):

    # prop masses
    Mp = rocket_equation_reversed(delta_v, isp, m_i)
    Mp_2 = rocket_equation_reversed(delta_v_2, isp, m_i)
    Mp_list.append(Mp)
    Mp_list_2.append(Mp_2)

    # times
    ttt = delta_v/(T[i]/m_i)
    # print("time:", ttt)
    mass_flow = T[i]/(isp*g0)
    t1t = Mp/mass_flow
    # print("time-t", t1t)

# print all masses for comparison
print(Mp_list)
print(Mp_list_2)

# volume calculation using EPSS 1CK (the best option)
v_epss = Mp_list[-1] * 1000 / rho[-1]


# Electric propulsion (non-impulsive manoeuvres)
# Check whether electric is possible or not
# ----------------------------
# options for electric, with their specific impulses and thrust levels
# PPTCUP, mu-CAT, PUC, Resistojet System, electrospray, MPACS, enpulsion, FFPT
isps_2 = [655, 3000, 70, 99, 800, 827, 3000, 1600, 2000]
T_2 = [40/(10**6), 50/(10**6), 5/1000, 100/1000, 0.7/1000, 80/(10**6),0.35/1000, 0.35/1000, 235/(10**6)]

# assuming circular to circular orbit (LEO to GEO) can also be used to check smaller changes
r = 42164
r0 = 200 + 6378

t_list = []
for i in range(len(isps_2)):
    t = Transfer_Time(m_i, isps_2[i], T_2[i] / 1000, r, r0)
    t /= 60
    t /= 60
    t /= 24
    t_list.append(t)

# print list with times for comparison
print("transfer times: ", t_list)
