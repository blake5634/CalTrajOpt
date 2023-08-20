#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import brl_data.brl_data as bd
import sys
import csv

fname = sys.argv[1]

print('opening ', fname)
fp = open(fname, newline='')

data = csv.reader(fp,delimiter=',', quotechar='"')

dt = 0.1
t = []
x = []
y = []
z = []
for row in data:
    print (row)
    i = int(row[0])
    x.append(float(row[1]))
    y.append(float(row[2]))
    z.append(float(row[3]))
    t.append(i*dt)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot3D(x,y,z, 'blue')
plt.show()
