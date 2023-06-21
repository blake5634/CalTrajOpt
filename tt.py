import c2to as cto
import matplotlib.pyplot as plt


cto.configure()  # read in config parameters

p1 = cto.point1D(0,0)
p2 = cto.point1D(3,1)


p1.x = -0.333
p1.v = -1

p2.x = -1.0
p2.v =  -0.333

tr = cto.trajectory1D(p1,p2)

#dt = 1.0

#tr.compute(dt)
tr.constrain_A()   # adjust dt for Amax
t,x,v,a = tr.timeEvolution()

plt.figure()
plt.plot(t,x)
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-3,3])
plt.title('Position v Time')
plt.grid(True)


plt.figure()
plt.plot(t,v)
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-5,5])
plt.title('Velocity')
plt.grid(True)

plt.figure()
plt.plot(t,a)
ax = plt.gca()
ax.set_xlim([0,3])
ax.set_ylim([-10,10])
plt.title('Acceleration')
plt.grid(True)


plt.show()
