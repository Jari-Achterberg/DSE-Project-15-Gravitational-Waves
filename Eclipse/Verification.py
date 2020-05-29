import numpy as np
import matplotlib.pyplot as plt
import Reference_frame

#%%
#Verification of the period formula: The result of this when using the Sun gravitational parameter and Earth distance should be around 365 days.
V_year = Reference_frame.T(Reference_frame.mu_sun,Reference_frame.a_earth)/3600/24

#%%
# The distance of the Earth relative to the sun in periapsis should be around 1.47*10^8 [km]
V_per_es = Reference_frame.r(Reference_frame.a_earth,Reference_frame.ecc_earth,0)

#%%


i = np.linspace(0,2*np.pi,200)
ecc = 1#Reference_frame.ecc_earth

M_c = Reference_frame.M_norm(i,ecc,0)
plt.plot(i,M_c)
plt.legend()
plt.show()

#%%


def rec_bis(N_max):
    '''
    :param a_0: Upper limit (has to be greater than zero)
    :param b_0: Lower limit (has to be smaller than zero)
    :param N_max: Number of iterations
    :return: Returns the value of the function using recursive bisection. In this case, E
    '''
    y = []
    a = 2*np.pi  # Upper limit
    b = -2*np.pi  # Lower limit
    for i in range(N_max):
        c = (a + b) / 2
        if Reference_frame.M_norm(a, Reference_frame.ecc_earth, np.pi-0.001) * Reference_frame.M_norm(c, Reference_frame.ecc_earth, np.pi-0.001) <= 0:
            b = c
        else:
            a = c
        y.append(c)
    y = np.array(y)
    return y

a = rec_bis(70)
print(a)
plt.plot(rec_bis(70))