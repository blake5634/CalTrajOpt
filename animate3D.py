#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D
import brl_data.brl_data as bd
import sys
#import csv

fname = sys.argv[1]

print('opening ', fname)
df = bd.datafile('', '','')
df.set_folders('','')
df.open('r',tname=fname)

dt = 0.1
t = []
x = []
y = []
z = []
for row in df.reader:
    print (row)
    i = int(row[0])
    x.append(float(row[1]))
    y.append(float(row[2]))
    z.append(float(row[3]))
    t.append(i*dt)


# References
# https://gist.github.com/neale/e32b1f16a43bfdc0608f45a504df5a84
# https://towardsdatascience.com/animations-with-matplotlib-d96375c5442c
# https://riptutorial.com/matplotlib/example/23558/basic-animation-with-funcanimation

# ANIMATION FUNCTION
def func(num, dataSet, line):
    # NOTE: there is no .set_data() for 3 dim data...
    line.set_data(dataSet[0:2, :num])
    line.set_3d_properties(dataSet[2, :num])
    return line

dataSet = np.array([x, y, z])
numDataPoints = len(z)

# GET SOME MATPLOTLIB OBJECTS
fig = plt.figure()
ax = Axes3D(fig)

# NOTE: Can't pass empty arrays into 3d version of plot()
line = plt.plot(dataSet[0], dataSet[1], dataSet[2], lw=2, c='g')[0] # For line plot

# AXES PROPERTIES]
# ax.set_xlim3d([limit0, limit1])
ax.set_xlabel('X(t)')
ax.set_ylabel('Y(t)')
ax.set_zlabel('Z(t)')
ax.set_title('3D phase space "optimal" trajectory')

# Creating the Animation object
line_ani = animation.FuncAnimation(fig, func, frames=numDataPoints, fargs=(dataSet,line), interval=25, blit=False)
#line_ani.save(r'AnimationNew.mp4')


plt.show()
