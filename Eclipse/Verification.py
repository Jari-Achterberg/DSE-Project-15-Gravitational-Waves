import numpy as np
import matplotlib.pyplot as plt
#import Reference_frame

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

#%%
t = np.linspace(1,10,10)
t_adj = t-5*np.trunc(t/5)
M = np.sqrt(8/10**3)*t_adj #E-e*sin(E)

j = np.array([1,2,3])>np.array([3,2,1])
k = j*np.array([1,2,3])
j_i = np.invert(j)

#%%
def pos_sun_earth(t,N_max): #Position of the Earth relative to the Sun after t seconds without taking only Earth orbit plane inclination into account
    '''
    :param t: Time after the pass through the pericentre for which the position of the Earth is obtained [s]
    :param N_max: Number of iterations for the recursive bisection
    :return: Cartesian coordinates of the position of the Earth relative to the Sun [km]
    '''
    t_adj = np.array(t)
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
        a = 2*np.pi*np.ones(len(t_adj))  # Upper limit
        b = -2*np.pi*np.ones(len(t_adj))  # Lower limit
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
   
    r_pos = r(a_earth,ecc_earth,theta_p) # Get the module of the r vector
    
    v_fin = np.array([r_pos*np.cos(theta_p)*np.cos(i_earth),r_pos*np.sin(theta_p),r_pos*np.cos(theta_p)*np.sin(i_earth)]) #Get the position in cartesian coordinates
    v_fin = -1*v_fin #Position of the sun relative to the Earth in a non-inclined and no-spinning axis of Earth 
    v_fin = np.matmul(T_polar,v_fin) #Calculation of the rotated reference frame, accounting for the polar inclination

    T_rot = np.array([[np.cos(Omega_t*t_adj),np.sin(Omega_t*t_adj),0],[-np.sin(Omega_t*t_adj),np.cos(Omega_t*t_adj),0],[0,0,1]]) #Transformation matrix to account for the rotation (days) of Earth
    v_fin = np.matmul(T_rot,v_fin)
    return v_fin
