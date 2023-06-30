import c2to as cto
import matplotlib.pyplot as plt


p1 = cto.point2D(0,0)
p2 = cto.point2D(3,1)


p1.x =  -2
p1.v = -1

p2.x =  1.5
p2.v =  +1.5

tr = cto.trajectory2D(p1,p2)

dt = 1.0

tr.compute(dt)
tr.constrain_A()
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
