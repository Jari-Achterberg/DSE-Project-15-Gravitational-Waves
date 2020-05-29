import numpy as np
from math import log10

#Constant
k = 1.38064852*10**-23  #[m^2*kg/(K*s^2)

#INPUTS REQUIRED
P   = 1    #[W]    output power of antenna (s/c) -->
L_l = 1    #[-]    Feed Line Loss  (s/c)
G_t = 1    #[-]    Antenna Gain (s/c)
e_t = 1    #[deg]  Pointing offset angle
a_half = 1 #[deg]  Antenna half-power beamwidth
d   = 1    #[m]    Distance s/c to ground station
L_a = 1    #[-]    Transmission path losses (from atmospheric attenuation and rain attenuation GRAPHS)
G_r = 1    #[-]    Antenna gain
A_r = 1    #[m^2]  Antenna effective area (ground)
D_r = 1    #[m]    Diameter antenna area (ground)
n_r = 1    #[-]    Antenna efficiency (ground)
L_r = 1    #[-]    Reception feeder loss (ground)
f   = 1    #[Hz]   Signal frequency
T_s = 1    #[K]    Noise Temperature
R   = 1    #[b/s]  Data Rate (FROM REQUIREMENTS)

#Calculated inputs
lam = 3*10**8/f             #[m]  Wavelength of signal
L_s = (lam/(4*np.pi*d))**2  #[-]  Space loss
L_pr = -12*(e_t/a_half)**2  #[-]  Antenna pointing loss

#Intermediate steps for sanity check
N_0 = k*T_s                          #[N*m]   Noise spectral density
EIRP = P*L_l*G_t                     #[W]
W_f  = EIRP/(4*np.pi*d**2)           #[W/m^2] Flux @ ground
C =  P*L_l*G_t*L_a*G_r*L_s*L_pr*L_r  #[W]     Power Received
E_b = C/R                            #[J/b]   Received Energy per bit

#Standard units
SNR = P * L_l * G_t * L_a * G_r * L_s * L_pr * L_r / (R * k * T_s)
print(SNR)

# #In Decibels
# SNR_d = P + L_l + G_t + L_a + G_r + L_s + L_pr + L_r + 228.6 - 10*log10(R) - 10*log10(T_s)
# print(SNR_d)