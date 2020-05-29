import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%%
ecc_earth = 0.01670861 #Earth eccentricity
a_earth = 149.6*10**6 # Semi-major axis of Earth [km]
mu_sun = 1.3271244*10**11 # Gravitational parameter of the sun [km^3s^-2]
i_earth = 7.155/180*np.pi # Earth inclination w.r.t the sun[rad]
i_polar = 23.43664/180*np.pi # Earth polar axis inclination relative to orbital plane
day = 86164.09053083288 # Sideral day duration

#%%
Or_p = 20000 #Number of points to modelate the orbit
N_max = 70 #Number of iterations for convergence of the recursive bisection

#%%

def T(mu,a): #orbital period
    '''
    :param mu: Gravitational parameter [km^3s^-2]
    :param a: Semi-major axis [km]
    :return: Orbital period [s]
    '''
    return 2*np.pi*np.sqrt(a**3/mu)

def r(a,e,theta):
    '''
    :param a: Semi-major axis [km]
    :param e: eccentricity
    :param theta: True anomaly [rad]
    :return: Module of the radius that goes from the centre of the planet/star [km]
    '''
    return a*(1-e**2)/(1+e*np.cos(theta))

ipos_earth = np.array([r(a_earth,ecc_earth,0),0]) #Initial position of the Earth at pericentre

def M_norm(E,e,M):
    '''
    :param E: Eccentric anomaly [rad]
    :param e: Eccentricity
    :param M: Mean anomaly [rad]
    :return: A theoretical value of zero for the equation that solves for recursive bisection
    '''
    return E-e*np.sin(E)-M

def pos_earth_sun(t,N_max): #Position of the Earth relative to the Sun after t seconds
    '''
    :param t: Time after the pass through the pericentre for which the position of the Earth is obtained
    :param N_max: Number of iterations for the recursive bisection
    :return: Cartesian coordinates of the position of the Earth relative to the Sun.
    '''
    T_pos = T(mu_sun,a_earth) #Period
    t_adj = t-T_pos*np.trunc(t/T_pos) #Adjusted time taking away the posibility of more than one revolution

    M = np.sqrt(mu_sun/a_earth**3)*t_adj #E-e*sin(E)

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
            if M_norm(a, ecc_earth, M) * M_norm(c, ecc_earth, M) <= 0: 
                #If I just put '<', there is an error at M=np.pi. This is something that should be looked into if there's time
                b = c
            else:
                a = c
            y.append(c)
        y = np.array(y)
        return y

    E_list = rec_bis(N_max) #Get the history of calculated E values [rad]
    #plt.plot(E_list)
    #plt.show()
    
    theta_p = np.arccos((1-ecc_earth**2)/ecc_earth/(1-ecc_earth*np.cos(E_list[-1]))-1/ecc_earth) #Get the true anomaly [rad]
    if t_adj > T_pos/2: #Correct for angles larger than pi
        theta_p = 2*np.pi-theta_p
   
    r_pos = r(a_earth,ecc_earth,theta_p) #Get the position in cartesian coordinates
    return np.array([r_pos*np.cos(theta_p)*np.cos(i_earth),r_pos*np.sin(theta_p),r_pos*np.cos(theta_p)*np.sin(i_earth)]) 


x_e_s = np.linspace(0,T(mu_sun,a_earth)/365,Or_p) #Generate points in time that account for one orbit
y_e_s = []
for i in x_e_s:
    y_e_s.append(pos_earth_sun(i,N_max))
y_e_s = np.transpose(np.array(y_e_s)) #Cartesian vector that represents the position of the Earth relative to the sun [km]
    
y_s_e = -1*y_e_s #Position of the sun relative to the Earth in a non-inclined and no-spinning axis of Earth 
T_polar = np.array([[np.cos(i_polar),0,-np.sin(i_polar)],[0,1,0],[np.sin(i_polar),0,np.cos(i_polar)]]) #Transformation matrix to account for the tilt of the polar axis of Earth

y_s_eT = np.transpose(y_s_e) #Transpose of the matrix with the positions of the sun relative to Earth
y_s_erT = [] #Calculation of the rotated reference frame, accounting for the polar inclination
for i in range(Or_p): 
    y_s_erT.append(np.matmul(T_polar,y_s_eT[i]))
y_s_erT = np.array(y_s_erT)
y_s_er = y_s_erT.T


Omega_t = 2*np.pi/day #Rotational rate of the Earth


y_s_errT = [] #Calculation of the rotated reference frame, accounting for the polar inclination and Earth's rotation
for i in range(Or_p):
    T_rot = np.array([[np.cos(Omega_t*x_e_s[i]),np.sin(Omega_t*x_e_s[i]),0],[-np.sin(Omega_t*x_e_s[i]),np.cos(Omega_t*x_e_s[i]),0],[0,0,1]]) #Transformation matrix to account for the rotation (days) of Earth
    y_s_errT.append(np.matmul(T_rot,y_s_erT[i]))
y_s_errT = np.array(y_s_errT)
y_s_err = y_s_errT.T


fig = plt.figure()
ax = fig.gca(projection='3d')
#ax.plot(y_s_e[0],y_s_e[1],y_s_e[2], label='Sun to Earth')
#ax.plot(y_s_er[0],y_s_er[1],y_s_er[2], label='Sun to Earth inclined')
ax.plot(y_s_err[0],y_s_err[1],y_s_err[2], label='Sun to Earth inclined and days')
ax.legend()
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
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
        if M_norm(a, ecc_earth, np.pi-0.001) * M_norm(c, ecc_earth, np.pi-0.001) <= 0:
            b = c
        else:
            a = c
        y.append(c)
    y = np.array(y)
    return y

a = rec_bis(70)
print(a)
plt.plot(rec_bis(70))