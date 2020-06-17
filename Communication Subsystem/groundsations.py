""" Very basic program for estimating Ground Coverage of
    satellites in relation to Ground Station positions """

import matplotlib.pyplot as plt
import numpy as np


# =================== GROUND STATION POSITIONS ===================

# Convert Lat-Long coords to Decimal Degrees
def dd(Deg, Min, Sec):
    return Deg + (Min / 60.) + (Sec / 3600.)

# Importing GS data from CSV file
filename = 'groundstations.csv'
gs = np.genfromtxt(filename, delimiter=',', dtype='str')

# Tables to store values for plotting
ddCoordsTab = []
gsNamesTab = []
gsColorsTab = ['royalblue', 'darkorange', 'yellowgreen', 'firebrick', 'plum']

# Calculating dd coords for each GS
for i in range(1, len(gs)):
    lat_i = gs[i][1]
    long_i = gs[i][2]

    # Stripping extra " from strings
    lat_i = lat_i.strip('"')
    long_i = long_i.strip('"')

    # Splitting coords strings
    lat_i = lat_i.split()
    long_i = long_i.split()

    # Calculating coords in Decimal Degrees
    if len(lat_i) == 2:
        # Latitude
        Lat_deg_i = lat_i[0]
        Lat_deg_i = float(Lat_deg_i.strip('°'))

        Lat_min_i = lat_i[1]
        Lat_min_i = float(Lat_min_i.strip("'"))

        Lat_sec_i = 0.

        # Longitude
        Long_deg_i = long_i[0]
        Long_deg_i = float(Long_deg_i.strip('°'))

        Long_min_i = long_i[1]
        Long_min_i = float(Long_min_i.strip("'"))

        Long_sec_i = 0.

    elif len(lat_i) == 3:
        # Latitude
        Lat_deg_i = lat_i[0]
        Lat_deg_i = float(Lat_deg_i.strip('°'))

        Lat_min_i = lat_i[1]
        Lat_min_i = float(Lat_min_i.strip("'"))

        Lat_sec_i = lat_i[2]
        Lat_sec_i = float(Lat_sec_i.strip("'"))

        # Longitude
        Long_deg_i = long_i[0]
        Long_deg_i = float(Long_deg_i.strip('°'))

        Long_min_i = long_i[1]
        Long_min_i = float(Long_min_i.strip("'"))

        Long_sec_i = long_i[2]
        Long_sec_i = float(Long_sec_i.strip("'"))

    Lat_i = dd(Lat_deg_i, Lat_min_i, Lat_sec_i)
    Long_i = dd(Long_deg_i, Long_min_i, Long_sec_i)

    # Recording in tables
    ddCoordsTab.append((Long_i, Lat_i))
    gsNamesTab.append(gs[i][0])

# ====================== SATELLITE POSITIONS ======================

# Parameters
satCoverage = 70.  # [deg] longitude
satDistance = 120.  # [deg] longitude
centerSat0 = -110.  # [deg] longitude


# Calculate boundaries of satellite coverage.
# Assumed that 75 [deg] coverage is favorable for comms at GEO.
def satBounds(center):
    left = float(center) - satCoverage
    right = float(center) + satCoverage
    return left, float(center), right


# Check for boundaries that are larger than 180 [deg] to plot the overlap.
def checkOverlap(bound):
    if abs(bound) > 180:
        overlap = float(abs(bound)) - 180.
    else:
        overlap = 0.
    return float(overlap)


# Position satellites at 120 [deg] from each other in longitude.
# The "center" satellite is taken as reference ("satellite 0").
# One is placed left ("satellite 1") and one is placed
# right of this center satellite ("satellite 2").
def geoSats(centerSat0):
    centerSat1 = float(centerSat0) + satDistance
    centerSat2 = float(centerSat0) - satDistance
    return float(centerSat0), centerSat1, centerSat2


# Calculating satellite boundaries for each satellite
satBounds0 = satBounds(geoSats(centerSat0)[0])
satBounds1 = satBounds(geoSats(centerSat0)[1])
satBounds2 = satBounds(geoSats(centerSat0)[2])


# ======== CALCULATE [GROUND STATION]-[SATELLITE] DISTANCE ========
# source: https://math.stackexchange.com/questions/301410/how-to-calculate-distance-from-a-given-latitude-and-longitude-on-the-earth-to-a

def GStoSAT(long_gs, lat_gs, long_sat):
    Re = 6371  # [km]
    Rs = 35786 + Re  # [km]
    lat = abs(lat_gs) * np.pi / 180.
    long = abs(long_gs - long_sat) * np.pi / 180.
    dist_squared = (Rs - Re * np.cos(long) * np.cos(lat)) ** 2 + (Re * np.sin(long) * np.cos(lat)) ** 2 + (
                Re * np.sin(lat)) ** 2
    dist = np.sqrt(dist_squared)
    return dist


# =========================== PLOTTING ===========================

# Basic plot parameters
fig = plt.figure(figsize=(7, 5))
fig.suptitle('Ground station VS Satellite positioning', fontsize=18, fontweight='bold')

axs = fig.add_subplot(1, 1, 1)

axs.set_ylabel('Longitude [deg]', fontsize=14)
axs.set_xlabel('Latitude [deg]', fontsize=14)
axs.tick_params(axis='both', which='major', labelsize=10)
# axs.grid(b = None, which = 'both', axis = 'both')

plt.xlim(-180, 180)
plt.ylim(-75, 75)
plt.axhline(0, color='black', lw='1')

# Plotting satellite 0
plt.axvline(satBounds0[0], color='deepskyblue', ls='dashed', lw=1)
plt.axvline(satBounds0[2], color='deepskyblue', ls='dashed', lw=1)
plt.axvspan(satBounds0[0], satBounds0[2], color='gray', alpha=0.2)
plt.plot(geoSats(centerSat0)[0], 0, color='deepskyblue', marker='o', markersize=13, label='Satellite A')

# Plotting satellite 1
plt.axvline(satBounds1[0], color='tomato', ls='dashed', lw=1)
plt.axvline(satBounds1[2], color='tomato', ls='dashed', lw=1)
plt.axvline((-180 + checkOverlap(satBounds1[2])), color='tomato', ls='dashed', lw=1)
plt.axvspan(satBounds1[0], satBounds1[2], color='gray', alpha=0.2)
plt.axvspan(-180, (-180 + checkOverlap(satBounds1[2])), color='gray', alpha=0.2)
if geoSats(centerSat0)[1] > 180:
    plt.plot((-180 + checkOverlap(geoSats(centerSat0)[1])), 0, color='tomato', marker='o', markersize=13,
             label='Satellite B')
else:
    plt.plot(geoSats(centerSat0)[1], 0, color='tomato', marker='o', markersize=13, label='Satellite B')

# Plotting satellite 2
plt.axvline(satBounds2[0], color='forestgreen', ls='dashed', lw=1)
plt.axvline(satBounds2[2], color='forestgreen', ls='dashed', lw=1)
plt.axvline((180 - checkOverlap(satBounds2[0])), color='forestgreen', ls='dashed', lw=1)
plt.axvspan(satBounds2[0], satBounds2[2], color='gray', alpha=0.2)
plt.axvspan((180 - checkOverlap(satBounds2[0])), 180, color='gray', alpha=0.2)
if geoSats(centerSat0)[2] < -180:
    plt.plot((180 - checkOverlap(geoSats(centerSat0)[2])), 0, color='forestgreen', marker='o', markersize=13,
             label='Satellite C')
else:
    plt.plot(geoSats(centerSat0)[2], 0, color='forestgreen', marker='o', markersize=13, label='Satellite C')

# Plotting ground stations
for i in range(len(gsNamesTab)):
    plt.plot(ddCoordsTab[i][0], ddCoordsTab[i][1], color=gsColorsTab[i], marker='^', markersize=15, label=gsNamesTab[i])

legend = fig.legend(bbox_to_anchor=(1.22, 0.885), framealpha=1, frameon=True, fontsize=10)  # prop={'size': 8})

frame = legend.get_frame()

frame.set_facecolor('0.98')
frame.set_edgecolor('black')

plt.show()

print('Distance between satellites: ', satDistance, ' [deg]')
print('Satellite coverage: ', satCoverage, ' [deg]')
print()
print('Satellites positioned at:')
print('Sat 0: ', satBounds0[1], ' [deg]')
print('Sat 1: ', satBounds1[1], ' [deg]')
print('Sat 2: ', satBounds2[1], ' [deg]')

# Distance from Sat A to Kourou [1]
dist01 = GStoSAT(ddCoordsTab[1][0], ddCoordsTab[1][1], geoSats(centerSat0)[0])

# Distance from Sat B to Kourou [1]
dist11 = GStoSAT(ddCoordsTab[1][0], ddCoordsTab[1][1], geoSats(centerSat0)[1])

# Distance from Sat B to Villafranca [2]
dist12 = GStoSAT(ddCoordsTab[2][0], ddCoordsTab[2][1], geoSats(centerSat0)[1])

# Distance from Sat B to Redu [3]
dist13 = GStoSAT(ddCoordsTab[3][0], ddCoordsTab[3][1], geoSats(centerSat0)[1])

# Distance from Sat B to Kiruna [4]
# dist14 = GStoSAT(ddCoordsTab[4][0],ddCoordsTab[4][1],geoSats(centerSat0)[1])

# Distance from Sat C to Perth [0]
dist21 = GStoSAT(ddCoordsTab[0][0], ddCoordsTab[0][1], geoSats(centerSat0)[2])

print()
print('A-Kourou =      ', dist01, ' km')
print('B-Kourou =      ', dist11, ' km')
print('B-Villafranca = ', dist12, ' km')
print('B-Redu =        ', dist13, ' km')
# print(dist14)
print('C-Perth =       ', dist21, ' km')

