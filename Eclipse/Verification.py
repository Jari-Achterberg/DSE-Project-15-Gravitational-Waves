import numpy as np
import matplotlib.pyplot as plt

i = np.linspace(0,2*np.pi,200)
ecc_earth = 0.999

def M(E,e):
    return E-e*np.sin(E)-np.pi

M_c = M(i,ecc_earth)
plt.plot(i,M_c)
plt.show()
