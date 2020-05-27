import numpy as np
import matplotlib.pyplot as plt

ecc_earth = 0.1
M = 0.9702068506

def f1(E,e,M):
    '''
    This formula calculates the position of a satellite at a particular point in time.
    :param E: Eccentric anomaly
    :param e: Eccentricity
    :param M: Mean anomaly
    :return: Returns a theoretical value of zero so the equation solver underneath can solve the equation
    '''
    return E-e*np.sin(E)-M

def f(a_0,b_0,N_max):
    '''
    This method uses the recursive bisection to solve the equation
    :param a_0: Upper limit (has to be greater than zero)
    :param b_0: Lower limit (has to be smaller than zero)
    :param N_max: Number of iterations
    :return: Returns the value of the function. In this case, E
    '''
    y = []
    a=a_0 #Upper limit
    b=b_0 #Lower limit
    for i in range(N_max):
        c = (a+b)/2
        if f1(a,ecc_earth,M)*f1(c,ecc_earth,M)<0:
            b=c
        else:
            a=c
        y.append(c)
    y=np.array(y)
    return y

theta = np.arccos((1-ecc_earth**2)/ecc_earth/(1-ecc_earth*np.cos(f(2*np.pi,-2*np.pi,50)))-1/ecc_earth)
print(theta)

final = f(2*np.pi,-2*np.pi,50)

plt.plot(final)
plt.show()
print(final[-1])