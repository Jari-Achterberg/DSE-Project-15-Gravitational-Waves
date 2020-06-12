import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#%%
ecc_earth = 0.01670861 #Earth eccentricity
a_earth = 149.6*10**6 # Semi-major axis of Earth [km]
mu_sun = 1.3271244*10**11 # Gravitational parameter of the sun [km^3s^-2]
i_earth = 0 /180*np.pi # Earth inclination w.r.t the sun[rad]
i_polar = (23.43664-i_earth)/180*np.pi # Earth polar axis inclination relative to orbital plane
day = 86164.09053083288 # Sideral day duration

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
T_polar = np.array([[np.cos(i_polar),0,-np.sin(i_polar)],[0,1,0],[np.sin(i_polar),0,np.cos(i_polar)]]) #Transformation matrix to account for the tilt of the polar axis of Earth
Omega_t = 2*np.pi/day #Rotational rate of the Earth

def M_norm(E,e,M):
    '''
    :param E: Eccentric anomaly [rad]
    :param e: Eccentricity
    :param M: Mean anomaly [rad]
    :return: A theoretical value of zero for the equation that solves for recursive bisection
    '''
    return E-e*np.sin(E)-M

#%%
Or_p = 10000
#Or_p = int(np.around(T(mu_sun,a_earth)/60))#Number of points to modelate the orbit
N = 35 #Number of iterations for convergence of the recursive bisection

#%%
def pos_sun_earth(t,N_max): #Position of the Earth relative to the Sun after t seconds without taking only Earth orbit plane inclination into account
    '''
    :param t: Time after the pass through the pericentre for which the position of the Earth is obtained [s] (Has to be an array)
    :param N_max: Number of iterations for the recursive bisection
    :return: Cartesian coordinates of the position of the Earth relative to the Sun [km]
    '''

    T_pos = T(mu_sun,a_earth) #Period
    N_p = np.trunc(t/T_pos) #Number of periods
    t_adj = t-T_pos*N_p #Adjusted time taking away the posibility of more than one revolution

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
            TrueFalse = M_norm(a, ecc_earth, M) * M_norm(c, ecc_earth, M) <= np.zeros(len(t_adj)) 
                #If I just put '<', there is an error at M=np.pi. This is something that should be looked into if there's time
            b = c*TrueFalse+b*np.invert(TrueFalse)
            a = a*TrueFalse+c*np.invert(TrueFalse)
            y.append(c)
        y = np.array(y)
        return y

    E_list = rec_bis(N_max) #Get the history of calculated E values [rad]
    #plt.plot(E_list)
    #plt.show()
    
    TrueFalse = t_adj > T_pos/2*np.ones(len(t_adj)) #Correct for angles larger than pi
    theta_p = np.arccos((1-ecc_earth**2)/ecc_earth/(1-ecc_earth*np.cos(E_list[-1]))-1/ecc_earth) #Get the true anomaly [rad]
    theta_p = TrueFalse*(2*np.pi-theta_p)+np.invert(TrueFalse)*theta_p
   
    r_pos = r(a_earth,ecc_earth,theta_p) # Get the module of the r vector
    
    v_fin = np.array([r_pos*np.cos(theta_p)*np.cos(i_earth),r_pos*np.sin(theta_p),r_pos*np.cos(theta_p)*np.sin(i_earth)]) #Get the position in cartesian coordinates
    v_fin = -1*v_fin #Position of the sun relative to the Earth in a non-inclined and no-spinning axis of Earth 
    v_fin = np.matmul(T_polar,v_fin) #Calculation of Sun position relative to Earth accounting for polar inclination

    T_rot = [] #Calculate the matrices that account for the Earth's spinning 
    for i in range(len(t_adj)):
        T_rot.append(np.array([[np.cos(Omega_t*t[i]),np.sin(Omega_t*t[i]),0],[-np.sin(Omega_t*t[i]),np.cos(Omega_t*t[i]),0],[0,0,1]])) #Transformation matrix to account for the rotation (days) of Earth
    T_rot = np.array(T_rot)
    
    v = []
    for i in range(len(t_adj)):
        c = np.matmul(T_rot[i],v_fin.T[i])
        v.append(c)
    v_fin = np.array(v)
    return v_fin.T

def time(date):
    '''
    Parameters
    ----------
    date : TYPE 2000-01-01T00:00:01
        DESCRIPTION. year-month-dayThour:miutes:seconds
    Returns
    -------
    TYPE float
        DESCRIPTION. Returns the amount of seconds that passed from the input and the 01-01-2000
    '''
    return (np.datetime64(date)-np.datetime64('2000-01-01')).astype("float")


def pos_sun_earth_time(begin,end,res,N_max): #Generate points in time that account for one orbit
    '''
    Parameters
    ----------
    begin : TYPE 2000-01-01T00:00:01
        DESCRIPTION. Begin date for the orbital calculations. 
    end : TYPE 2001-01-01T00:00:01
        DESCRIPTION. End date for the orbital calculations. In this case it would be a year orbit
    res : TYPE Int
        DESCRIPTION. Represents the size of the output vector 
    N_max : TYPE Int
        DESCRIPTION. Represents the number of iterations for the recursive bisection     
    Returns
    -------
    TYPE array
        DESCRIPTION. Returns a vector for the x, y and z points of the sun with respect to the Earth for every moment in time specified.    
    '''
    b = time(begin)
    e = time(end)
    x_s = np.linspace(b,e,res)
    f = pos_sun_earth(x_s,N_max)
    
    return np.array([f[0],f[1],f[2],x_s])


#%%


def pos_sun_sat(begin,end,res,N_max):
    '''
    Parameters
    ----------
    begin : TYPE 2000-01-01T00:00:01
        DESCRIPTION. Begin date for the orbital calculations. 
    end : TYPE 2001-01-01T00:00:01
        DESCRIPTION. End date for the orbital calculations. In this case it would be a year orbit
    res : TYPE Int
        DESCRIPTION. Represents the size of the output vector 
    N_max : TYPE Int
        DESCRIPTION. Represents the number of iterations for the recursive bisection     
    Returns
    -------
    TYPE array
        DESCRIPTION. Returns a 5*res array for the x, y and z points of the Sun with respect to the Satellite, the 
        vertical angle of the satellite wrt the sun and the angle of the Sun in the x-z plane for every moment in 
        time specified.    
    '''
    f = pos_sun_earth_time(begin,end,res,N_max)
    vec = np.matmul(np.array([[0,1,0],[0,0,1],[1,0,0]]),f[0:3,:])
    angle_sun = np.array([np.arctan(vec[1]/np.sqrt(vec[0]**2+vec[2]**2)),np.arctan2(vec[0],vec[2])])
    return np.array([vec[0],vec[1],vec[2],angle_sun[0],angle_sun[1],f[3]])
    
#%%

#y_sat = pos_sun_sat('2000-01-01T00:00:00','2001-01-01T00:00:00',Or_p,N) #Function for Sun position relative to satellite
#y_earth = pos_sun_earth_time('2000-01-01T00:00:00','2001-01-01T00:00:00',Or_p,N) #Function for Sun position relative to Earth
'''
y_earthT = y_earth.T
f = open("output.dat","w")
for i in range(len(y_earthT)):
    line = str(y_earthT[i,0])+","+str(y_earthT[i,1])+","+str(y_earthT[i,2])+","+str(y_earthT[i,3])+"\n"
    f.write(line)
f.close()
'''
y_earth = pos_sun_earth_time('2000-01-01T00:00:00','2001-01-01T00:00:00',Or_p,N)
'''
y_s = np.matmul(np.array([[0,1,0],[0,0,1],[1,0,0]]),y_earth)
angle_sun_sat = np.array([np.arctan(y_s[1]/np.sqrt(y_s[0]**2+y_s[2]**2)),np.arctan2(y_s[0],y_s[2])])
'''
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(y_earth[0],y_earth[1],y_earth[2], label='Sun to Earth inclined and days')
ax.legend()
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
plt.show()