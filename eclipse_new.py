import numpy as np

d = 149597870  # km, average distance from Earth to Sun
r_sun = 696340  # km, Sun radius
r_earth = 6378  # km, Earth radius
e = 0  # GEO eccentricity
i = 0  # rad, inclination of the orbit
omega = 0  # rad, right ascension of the ascending node
om = 0  # rad, argument of periapsis

a = 42164  # km, semi-major axis
p = 42164  # km, semi-minor axis


def sat_2d_pos(theta):
    """
    This function returns position of the satellite
    :param theta: true anomaly
    :return: position coordinates
    """
    r_sat = a * (1 - e**2) / (1 + e * np.cos(theta))
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

def vector_magnitude(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

### Transformation parameters

l1 = np.cos(omega) * np.cos(om) - np.sin(omega) * np.sin(om) * np.cos(i)
l2 = -np.cos(omega) * np.sin(om) - np.sin(omega) * np.cos(om) * np.cos(i)
m1 = np.sin(omega) * np.cos(om) + np.cos(omega) * np.sin(om) * np.cos(i)
m2 = -np.sin(omega) * np.sin(om) + np.cos(omega) * np.cos(om) * np.cos(i)
n1 = np.sin(om) * np.sin(i)
n2 = np.cos(om) * np.sin(i)

transformation_parameter = np.array([[l1, l2],
                                    [m1, m2],
                                    [n1, n2]])

# sat_2d_pos = sat_position(0.8421)
# sat_3d_pos = np.dot(transformation_parameter, xi_eta(sat_2d_pos))  # km, satellite position in 3D


def sat_3d_position(sat_2d_position):
    """
    Satellite position in 3D Cartesian
    :param sat_2d_position: satellite position in 2D Cartesian
    :return: position in 3D
    """
    return np.dot(transformation_parameter, xi_eta(sat_2d_position))

sun_r = np.array([[d],  # Sun position vector
                  [0],
                  [0]], dtype='int64')


def big_psi(sun_pos, sat_3d_pos):
    """
    Angle between vector to Sun and vector to satellite
    :param sun_pos: Sun position vector
    :param sat_3d_pos: satellite position in 3D Cartesian
    :return: angle in radians
    """
    return np.arccos(np.dot(sun_pos.T, sat_3d_pos) / (vector_magnitude(sun_pos[0], sun_pos[1], sun_pos[2]) * vector_magnitude(sat_3d_pos[0], sat_3d_pos[1], sat_3d_pos[2])))


def a_sat(sat_3d_pos, psi):
    """
    Component of satellite vector, perpendicular to direction towards the Sun
    :param sat_3d_pos: satellite position in 3D Cartesian
    :param psi: angle between vector to Sun and vector to satellite
    :return: perpendicular component
    """
    return vector_magnitude(sat_3d_pos[0], sat_3d_pos[1], sat_3d_pos[2]) * np.sin(psi)


### Long-term eclipse
def unit_sun_r(sun_pos):
    """
    Sun position unit vector
    :param sun_r: Sun position vector
    :return: unit vector
    """
    return sun_pos / vector_magnitude(sun_pos[0], sun_pos[1], sun_pos[2])

n_vector = np.array([[np.sin(i) * np.sin(omega)],
                     [-np.sin(i) * np.cos(omega)],
                     [np.cos(i)]])


def check_if_in_shadow(psi, a_sat_vector, sun_pos):
    """
    This function checks if the satellite is in shadow cone
    :param psi: angle between direction to Sun and direction to satellite
    :param a_sat_vector: component of satellite vector, perpendicular to direction towards the Sun
    :param sun_pos: sun position vector
    :return:
    """
    dot_prod = np.dot(a * n_vector.T, unit_sun_r(sun_pos))

    check = np.zeros((len(a_sat_vector)))
    if np.cos(psi) < 0 and a_sat_vector < r_earth and dot_prod <= r_earth:
        check = True

    return check


theta_range = np.linspace(0, 2 * np.pi, 100, endpoint=True)

two_d_pos_tuple = sat_2d_pos(theta_range)

two_d_pos_range = np.vstack((two_d_pos_tuple[0], two_d_pos_tuple[1]))

three_d_pos_range = sat_3d_position(two_d_pos_range)

angle = big_psi(sun_r, three_d_pos_range)

perpendicular_comp = a_sat(three_d_pos_range, angle)

print(perpendicular_comp < r_earth)
# print(np.dot(a*n_vector.T, unit_sun_r(sun_r)))
print(n_vector)