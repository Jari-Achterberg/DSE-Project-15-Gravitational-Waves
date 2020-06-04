import numpy as np
import matplotlib.pyplot as plt


J = 1367 # W/m^2
c = 3e8 # m/s
q = 0.88
height = 0.30
width = 0.20
length = 0.20


#angle between the sun and the satellite in a specific reference system
alphas = np.arange(0, 2 * np.pi, 0.001 * np.pi)#sun around axis


betas1 = np.arange(0.122, np.pi/6 , 0.001 * np.pi)
betas2 = np.arange(0.122, np.pi/6 , 0.001 * np.pi)[::-1]


betas = np.concatenate((betas1,betas2,betas1,betas2,betas1,betas2,betas1,betas2,betas1)) #sun around y-axis



def SolarArea(height,width,length,alpha,beta):
    v_width = width * np.cos(alpha)
    v_height =height * np.cos(beta)
    v_area = v_width * v_height

    z_width = length * np.sin(alpha)
    z_height = height * np.cos(beta)
    z_area = z_width * z_height

    o_width = length
    o_height = width
    o_area = o_width * o_height * np.sin(beta)

    total_area = v_area + z_area + o_area

    return total_area

Force_sun = []
for i in range(1152):
    A = SolarArea(height,width,length,alphas[i],betas[i])
    F = (J/c)*A*np.cos(np.pi/6)*(1+q)
    Force_sun.append(A)

plt.plot(alphas[:1152],Force_sun)
