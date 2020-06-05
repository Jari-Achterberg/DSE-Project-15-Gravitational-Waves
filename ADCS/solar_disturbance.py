import numpy as np
import math as m
import matplotlib.pyplot as plt


J = 1367 # W/m^2
c = 3e8 # m/s
q = 0.88
height = 0.30 #y-axis [m]
width = 0.20 #x-axis [m]
length = 0.20 #z-axis [m]
cg = 0.1,0.15,0.1 #assumed

#angle between the sun and the satellite in a specific reference system
alphas = np.arange(0, 2 * np.pi, 0.001 * np.pi)#sun around z-axis


betas1 = np.arange(0.122, np.pi/6 , 0.001 * np.pi)
betas2 = np.arange(0.122, np.pi/6 , 0.001 * np.pi)[::-1]


betas = np.concatenate((betas1,betas2,betas1,betas2,betas1,betas2,betas1,betas2,betas1,betas2)) #sun around y-axis
#betas = np.zeros(len(alphas))

def transformation(a):
    A = np.array([(0, 0, 1), (1, 0, 0)])
    B = np.array([a[0], a[1], a[2]])
    C = np.dot(A,B)
    return C

def SolarArea(height,width,length,alpha,beta):
    v_width = abs(width * np.cos(alpha)) #side facing earth
    v_height = abs(height * np.cos(beta))
    v_area = v_width * v_height
    v_cp = v_height/2, v_width/2

    z_length = abs(length * np.sin(alpha)) #side 90deg from earth perpendicular to orbit
    z_height = abs(height * np.cos(beta))
    z_area = z_length * z_height
    z_cp = -z_length/2 , -z_height/2

    o_length = length #side on "top" parallel to orbit
    o_width = width
    o_area = o_width * o_length * np.sin(beta)
    o_cp =o_length/2,  -o_width/2

    total_area = v_area + z_area + o_area
    cp = (v_cp[0]*v_area + z_cp[0]*z_area + o_cp[0]*o_area)/total_area, (v_cp[1]*v_area + z_cp[1]*z_area + o_cp[1]*o_area)/total_area
    cp = np.array(cp)
    return total_area, cp



Torque_sun = []
cg_2d = transformation(cg)
for i in range(len(betas)):
    A, cp = SolarArea(height,width,length,alphas[i],betas[i])
    F = (J/c)*A*(1+q)
    r = m.sqrt((cg_2d[0]- cp[0])**2+(cg_2d[1]-cp[1])**2)
    T = F*r
    Torque_sun.append(T)


plt.plot(alphas[:len(betas)],Torque_sun)


