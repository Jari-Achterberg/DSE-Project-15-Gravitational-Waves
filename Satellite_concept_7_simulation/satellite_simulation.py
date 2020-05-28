import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import math
from points_distance import closestDistanceBetweenLines as cdbl

###RUN THIS FOR PLOT IN SEPARATE WINDOW
#            %matplotlib qt


plt.close()
# attaching 3d axis to the figure
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# setting the axes properties
ax.set_xlim3d([-42164000, 42164000])
ax.set_xlabel('x')

ax.set_ylim3d([-42164000, 42164000])
ax.set_ylabel('y')

ax.set_zlim3d([-42164000, 42164000])
ax.set_zlabel('z')
ax.set_title('satellite constellation concept 5')


def deg_to_rad(rad):
    return math.pi * (rad / 180.)


class Orbit:
    SMA = i = lan = ta = 0

    def __init__(self, i=0., lan=0., aop=0., ta=0.):
        # Inclination [rad]
        self.i = i

        # Longitude of Ascending Node [rad]
        self.lan = lan

        # True Anomaly [rad]
        self.ta = ta

        # period
        self.period = (24 * 60 * 60)

    def getRadius(self, theta):
        r = 42164000.
        return r

    def getTheta(self, t):
        theta = self.ta + (2 * math.pi * (t / self.period))
        return theta

    def getPosition(self, t):
        theta = self.getTheta(t)
        r = self.getRadius(theta)

        current_inclination = np.sin(theta - self.lan) * np.sin(self.i)
        z = current_inclination * r

        square_radius = np.square(r) - np.square(z)

        radius_altered_by_inclination = np.sqrt(square_radius)

        x = np.cos(theta) * radius_altered_by_inclination
        y = np.sin(theta) * radius_altered_by_inclination
        return [x, y, z]


# first triangle
satellite1 = Orbit(i=deg_to_rad(0.))  # sending lasers
satellite2 = Orbit(i=deg_to_rad(0.), ta=deg_to_rad(120.))
satellite3 = Orbit(i=deg_to_rad(0.), ta=deg_to_rad(240.))

# second triangle
satellite4 = Orbit(i=deg_to_rad(1.), ta=deg_to_rad(180.))  # sending lasers
satellite5 = Orbit(i=deg_to_rad(1.), ta=deg_to_rad(60.))
satellite6 = Orbit(i=deg_to_rad(1.), ta=deg_to_rad(300.))

circle2 = plt.Circle((5, 5), 0.5, color='b', fill=False)


def draw_dot(pos):
    plt.plot([pos[0]], [pos[1]], [pos[2]], 'ro', markersize=20)


def draw_line(pos1, pos2):
    return plt.plot([pos1[0], pos2[0]], [pos1[1], pos2[1]], [pos1[2], pos2[2]], 'y')


def draw_figure(t):
    for i in range(len(ax.lines)-1, -1, -1):
        ax.lines[i].remove()

    sats = [satellite1, satellite2, satellite3, satellite4, satellite5, satellite6]
    positions = [sat.getPosition(t) for sat in sats]

    for pos in positions:
        draw_dot(pos)

    plt.plot([0], [0], [0], 'bo', markersize=50)

    draw_line(positions[0], positions[1])
    draw_line(positions[0], positions[2])
    draw_line(positions[1], positions[2])
    draw_line(positions[4], positions[5])
    draw_line(positions[3], positions[4])
    draw_line(positions[3], positions[5])


# Creating the Animation object
line_ani = animation.FuncAnimation(fig, draw_figure, np.linspace(0, 3 * 99900, 2000), interval=20, blit=False,
                                   repeat=True)

plt.show()

# THIS DECPRECATED CODE WAS USED FOR CALCULATING SATELLITE LASER INTERFERENCE
'''
t_data = np.arange(0,86400,10)
i_data = []
sats = [satellite1, satellite2, satellite3, satellite4, satellite5, satellite6]
for t in t_data:
    positions = [np.array(sat.getPosition(t)) for sat in sats]
    lines1 = []
    lines2 = []
    lines1.append([positions[0], positions[1]])
    lines1.append([positions[0], positions[2]])
    lines1.append([positions[1], positions[2]])
    lines2.append([positions[4], positions[5]])
    lines2.append([positions[3], positions[4]])
    lines2.append([positions[3], positions[5]])

    for line1 in lines1:
        for line2 in lines2:        
            a,b,r = cdbl(line1[0],line1[1],line2[0],line2[1])
            if r < 4000*(2/3):
                i_data.append(t)
                break
        break

fraction_of_intersection = (len(i_data)/len(t_data))

print(fraction_of_intersection)
'''


