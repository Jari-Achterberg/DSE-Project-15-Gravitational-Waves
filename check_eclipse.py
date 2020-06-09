import numpy as np
import matplotlib.pyplot as plt
import datetime as dt

d = 149597870  # km, average distance from Earth to Sun
r_sun = 696340  # km, Sun radius
r_earth = 6378  # km, Earth radius
earth_mass = 5.9722 * 10**24  # kg, Earth mass
grav_const = 6.67408 * 10**(-11)  # m3 * kg^-1 * s^2, gravitational constant
i = 0/180 * np.pi  # rad, inclination of the orbit
omega = 0  # rad, right ascension of the ascending node
om = 0  # rad, argument of periapsis
earth_day = 86164.09903691  # s, Earth rotation period

apogee = 42164  # km, apogee distance
perigee = 42164  # km, perigee distance
a_earth = (apogee + perigee) / 2  # km, semi-major axis
e_earth = (apogee - perigee) / (apogee + perigee)  # eccentricity
mu_earth = earth_mass * grav_const * 10**(-9)  # km3 * s2, Earth gravitational parameter


def sat_2d_pos(theta):
    """
    This function returns position of the satellite
    :param theta: true anomaly
    :return: position coordinates
    """
    r_sat = a_earth * (1 - e_earth**2) / (1 + e_earth * np.cos(theta))
    return r_sat, theta


def xi_eta(sat_pos):
    """
    This function returns satellite position in Cartesian coordinates in 2D
    :param sat_pos: satellite position in polar coordinates
    :return: satellite position in Cartesian coordinates
    """
    xi = sat_pos[0] * np.cos(sat_pos[1])
    eta = sat_pos[0] * np.sin(sat_pos[1])
    return np.array([xi,
                     eta])

l1 = np.cos(omega) * np.cos(om) - np.sin(omega) * np.sin(om) * np.cos(i)
l2 = -np.cos(omega) * np.sin(om) - np.sin(omega) * np.cos(om) * np.cos(i)
m1 = np.sin(omega) * np.cos(om) + np.cos(omega) * np.sin(om) * np.cos(i)
m2 = -np.sin(omega) * np.sin(om) + np.cos(omega) * np.cos(om) * np.cos(i)
n1 = np.sin(om) * np.sin(i)
n2 = np.cos(om) * np.sin(i)

transformation_parameter = np.array([[l1, l2],
                                    [m1, m2],
                                    [n1, n2]])


def sat_3d_position(sat_2d_position):
    """
    Satellite position in 3D Cartesian
    :param sat_2d_position: satellite position in 2D Cartesian
    :return: position in 3D
    """
    return np.matmul(transformation_parameter, sat_2d_position)


def vector_magnitude(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)


def big_psi(sun_pos, sat_3d_pos):
    """
    Angle between vector to Sun and vector to satellite
    :param sun_pos: Sun position vector
    :param sat_3d_pos: satellite position in 3D Cartesian
    :return: angle in radians
    """
    return np.arccos(np.clip(np.dot(sun_pos.T, sat_3d_pos) / (vector_magnitude(sun_pos) * vector_magnitude(sat_3d_pos)), -1, 1))


def a_sat(sat_3d_pos, psi):
    """
    Component of satellite vector, perpendicular to direction towards the Sun
    :param sat_3d_pos: satellite position in 3D Cartesian
    :param psi: angle between vector to Sun and vector to satellite
    :return: perpendicular component
    """
    return vector_magnitude(sat_3d_pos) * np.sin(psi)


### Long-term eclipse
def unit_vector(vector):
    """
    Sun position unit vector
    :param sun_r: Sun position vector
    :return: unit vector
    """
    return vector / vector_magnitude(vector)


n_vector = np.array([[np.sin(i) * np.sin(omega)],
                     [-np.sin(i) * np.cos(omega)],
                     [np.cos(i)]])


def check_if_in_shadow(psi, a_sat_vector, sun_pos):
    dot_prod = a_earth * np.dot(n_vector.T, unit_vector(sun_pos))
    check = False

    if np.cos(psi) < 0 and a_sat_vector < r_earth and dot_prod <= r_earth:
        check = True

    return check


sun_vecs = np.load("sun_vectors.npy")

sat_pos_polar = sat_2d_pos(0.5)
sat_2d_pos_cart = xi_eta(sat_pos_polar)
sat_3d_pos = sat_3d_position(sat_2d_pos_cart)

eclipse = np.zeros(len(sun_vecs))
count = 0

for sun_vec in sun_vecs:
    sun_vec = sun_vec[0: 3]
    angle = big_psi(sun_vec, sat_3d_pos)
    a_satellite = a_sat(sat_3d_pos, angle)
    check = check_if_in_shadow(angle, a_satellite, sun_vec)
    eclipse[count] = check
    count += 1


print("Done with iterations")

plt.plot(eclipse)
plt.show()