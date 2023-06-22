import c2to as cto
import matplotlib.pyplot as plt
import numpy as np


cto.configure()  # read in config parameters

    #def __init__(self,ix,iy,iz,ixd,iyd,izd):

#p1 = cto.point3D(8,7,6,4,3,2)
#p2 = cto.point3D(2,3,4,4,5,6)
p1 = cto.point3D(1,0,2,2,0,1)
p2 = cto.point3D(0,2,1,0,2,0)

tr = cto.trajectory3D(p1,p2)

#dt = 1.0

#tr.compute(dt)
tr.constrain_A()   # adjust dt for Amax
t,x,v,a = tr.timeEvolution()

#print(tmp)
#print('X shape:',tmp.shape)
#print('t[2].shape',tmp[2].shape, tmp[2])
#print('t[:][0].shape',tmp[:][0].shape)
#print('t[0][:].shape',tmp[0][:].shape)
#print('t[:][2].shape',tmp[:][2].shape, tmp[:][2])


# plot x,y,z trajectory/tracectories

axp = 2  #
axnames = ['X', 'Y', 'Z']

name = axnames[axp]

print('x shape (tt3D): ',x.shape)
plt.figure()
plt.plot(t,x[axp])
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-3,3])
plt.title('Position v Time: ' + axnames[axp])
plt.grid(True)


plt.figure()
plt.plot(t,v[axp])
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-5,5])
plt.title('Velocity vs Time: ' + axnames[axp])
plt.grid(True)

plt.figure()
plt.plot(t,a[axp])
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-(cto.AMAX+1),cto.AMAX+1])
plt.title('Acceleration: ' + axnames[axp])
plt.grid(True)


plt.show()
