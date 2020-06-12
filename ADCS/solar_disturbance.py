import scipy as np
import math as m
import matplotlib.pyplot as plt
from Eclipse.Reference_frame import pos_sun_sat
import time

# Constants
J = 1367 # W/m^2
c = 3e8 # m/s
q = 0.88
height = 0.30 #y-axis [m]
width = 0.20 #x-axis [m]
length = 0.20 #z-axis [m]
cg = 0.1,0.15,0.1 #assumed


# for one day
'''
Or_p = 86400
alphas = np.linspace(0, 2 * np.pi, num = 86400)#sun around z-axis
betas = 23.3*np.pi/180
t = np.arange(0, 86400)
'''
#for one year

Or_p = 525981  #Number of points to modelate the orbit
N = 35 #Number of iterations for convergence of the recursive bisection
y_sat = pos_sun_sat('2000-01-01T00:00:00','2001-01-01T00:00:00', Or_p, N)
alphas =  y_sat[4]
betas = y_sat[3]
t = y_sat[5]*(1/3600)*(1/24)

#Definitions

def transformation(a): #to project center of gravity in 2D
    A = np.array([(0, 0, 1), (1, 0, 0)])
    B = np.array([a[0], a[1], a[2]])
    C = np.dot(A,B)
    return C


def SolarArea(height,width,length,alpha,beta): #Computes which faces of the satellite are facing the sun and how mu radiation they get
    v_width = abs(width * np.cos(alpha)) #side facing earth
    v_height = abs(height * np.cos(beta))
    if alpha > np.pi/2 and alpha < 3/2*np.pi:
        v_area = 0
    else:
        v_area = v_width * v_height
    v_cp = v_height/2, v_width/2

    z_length = abs(length * np.sin(alpha))  # side 90deg from earth perpendicular to orbit
    z_height = abs(height * np.cos(beta))
    z_area = z_length * z_height
    z_cp = -z_length/2 , -z_height/2

    o_length = length   # side on "top" parallel to orbit
    o_width = width
    o_area = 0      # o_width * o_length * np.sin(beta)
    o_cp = 0, 0     # o_length/2,  -o_width/2

    total_area = v_area + z_area + o_area
    cp = (v_cp[0]*v_area + z_cp[0]*z_area + o_cp[0]*o_area)/total_area, (v_cp[1]*v_area + z_cp[1]*z_area + o_cp[1]*o_area)/total_area
    cp = np.array(cp)
    return total_area, cp

# Code

eclipse = np.load('/Users/jaria/PycharmProjects/DSE-Project-15-Gravitational-Waves/ADCS/eclipse_check.npy')

Torque_sun = []
Area_exposed = []
cg_2d = transformation(cg)
mz = 0
for i in range(len(alphas)):
    if eclipse[i] == 0:
        A, cp = SolarArea(height,width,length,alphas[i],betas[i])
        F = (J/c)*A*(1+q)
        r = m.sqrt((cg_2d[0]- cp[0])**2+(cg_2d[1]-cp[1])**2)
        T = F*r
        Torque_sun.append(T)
        Area_exposed.append(A)
    else:
        mz += 1
        A = 0
        T = 0
        Torque_sun.append(T)
        Area_exposed.append(A)
'''
time1 = []
for i in range(len(t)):
    t1 = t[i] + 1925383546
    t2 = time.ctime(t1)
    time1.append(t2)
'''
H = np.trapz(Torque_sun,t) #Angular Momentum
#figure

fig = plt.figure(figsize = (7,5))
fig.suptitle('Area', fontsize=18, fontweight='bold')

axs = fig.add_subplot(1,1,1)
axs.plot(t, Area_exposed, color = '#00A6D6', label = 'label')

print(sum(Area_exposed)/(len(Area_exposed)-mz))

axs.set_ylabel('Area [m^2]', fontsize=14)
axs.set_xlabel('Time [Days]', fontsize = 14)
axs.tick_params(axis='both', which='major', labelsize=10)
axs.grid(b = None, which = 'both', axis = 'both')
plt.show()

#In case legend is needed
'''
legend = fig.legend(loc= 0 , framealpha=1, frameon=True, fontsize=16)

fig.align_labels()

plt.xlabel('Time [s]', fontsize=16)
frame = legend.get_frame()

frame.set_facecolor('0.90')
frame.set_edgecolor('black')
'''
