#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt

def error(msg):
    print('Error: ')
    print(msg)
    quit()

N = 4
M = N*N  # number of distances

AMAX = 2  #  normalize for now.
DT_TEST = 2.0


NPC = 20  # number plot points per curve

costtype = 'time'
costtype = 'energy'

startrow = 3
startcol = 1

def configure(fp=None):
    global costtype, AMAX, DT_TEST, N,M, NPC, startrow, startcol
    if not fp:
        f = open('ctoConfig.txt','r')
    else:
        f = fp

    for line in f:
        parname, *rest = line.split(' ')
        parname = parname.strip()
        v = ' '.join(rest).strip()

        if parname == '#':
            next

        if parname =='costtype':
            if v not in ['time','energy']:
                error('path.configure: unknown cost type: '+v)
            costtype = v
        if parname == 'amax':
            AMAX = float(v)
        if parname == 'dt_test':
            DT_TEST = dt_test
        if parname == 'N':
            N = int(v)
            M = N*N
        if parname == 'NPC':
            NPC = int(v)
        if parname == 'startrow':
            startrow = int(v)
        if parname == 'startcol':
            startcol = int(v)




class grid1D:
    def __init__(self, N):
        self.N = N
        self.gr = []
        for i in range(self.N):
            row = []
            for j in range(self.N):
                row.append(point1D(i,j))
            self.gr.append(row)

    def __repr__(self):
        txt = ''
        for i in range(self.N):
            txt += '       '
            for j in range(self.N):
                txt +=  ' ({:5.2f},{:5.2f})'.format(self.gr[i][j].x, self.gr[i][j].v)
            txt += '\n'
        return txt

class point1D:
    def __init__(self):
        self.row = 0
        self.col = 0
        self.__init(self,self.row,self.col)

    def __init__(self,i,j):
        self.row = i
        self.col = j
        self.x =     2*j/(N-1) - 1
        self.v = -1*(2*i/(N-1) - 1)

    def __eq__(self,x):
        return self.row == x.row and self.col == x.col

    def __repr__(self):
        return ' ({:4.1f}, {:4.1f})'.format(self.x, self.v)
        #return ' ({:}, {:})'.format(self.row, self.col)


class point3D:

    def __init__(self,ix,iy,iz,ixd,iyd,izd):
        self.ix = ix
        self.iy = iy
        self.iz = iz
        self.ixd = ixd
        self.iyd = iyd
        self.izd = izd
        self.x =     2*ix/(N-1) - 1
        self.y =     2*iy/(N-1) - 1
        self.z =     2*iz/(N-1) - 1
        self.xd = -1*(2*ixd/(N-1) - 1)
        self.yd = -1*(2*iyd/(N-1) - 1)
        self.zd = -1*(2*izd/(N-1) - 1)
        self.ivect = [self.ix, self.iy, self.iz, self.ixd, self.iyd, self.izd]
        self.xvect = [self.x,self.y,self.z,self.xd,self.yd,self.zd]

    def __eq__(self,x):
        for i,s in self.ivect:
            if s != x.ivect[i]:
                return False
        return True

    def __repr__(self):
        st = ' ('
        for s in self.ivect:
            st += '{:3d}, '.format(s)
        st += ')'
        return st



class trajectory3D:
    def __init__(self,p1,p2):
        if str(type(p1)) != "<class 'c2to.point3D'>" or str(type(p2)) != "<class 'c2to.point3D'>":
            error('trajectory3D called with 1D point!')
        print('trajectory1D: p1,p2: ',p1,p2)
        self.p1 = p1
        self.p2 = p2
        self.a0 = [0,0,0]
        self.a1 = [0,0,0]
        self.a2 = [0,0,0]
        self.a3 = [0,0,0]
        self.computed = False
        self.constrained = False
        self.valid = True

    def compute(self,dt):
        p1 = self.p1
        p2 = self.p2
        for i in range(3):
            self.a0[i] = p1.xvect[i]
            self.a1[i] = p1.xvect[i+3]
        #define some constants
        b0 = dt
        b1 = dt*dt
        b2 = dt*dt*dt
        b3 = 2.0*dt
        b4 = 3.0*dt*dt
        for i in range(3):
            dx = p2.xvect[i]-p1.xvect[i]
            dv = p2.xvect[i+3]-p1.xvect[i+3]
            # a3 solution:
            self.a3[i] = (dv - (b3/b1)*(dx-p1.xvect[i+3]*b0))/(b4-b2*b3/b1)
            self.a2[i] = (dx - p1.xvect[i+3]*b0 - self.a3[i]*b2)/b1
        self.computed = True

    def x(self,t):
        rv = [0,0,0]
        for i in range(3):
            rv[i] = self.a0[i] + self.a1[i]*t +    self.a2[i]*t*t +     self.a3[i]*t*t*t
        return rv
    def xd(self,t):
        rv = [0,0,0]
        for i in range(3):
            rv[i] =     self.a1[i]   + 2.0*self.a2[i]*t  + 3.0*self.a3[i]*t*t
        return rv
    def xdd(self,t):
        rv = [0,0,0]
        for i in range(3):
            rv[i] =       2.0*self.a2[i]    + 6.0*self.a3[i]*t
        return rv


    def get_Amax(self,dt):
        if not self.computed:
            error('cant get_Amax() unless tr.computed is True')
        amax = -1000
        acc1 = self.xdd(0)
        acc2 = self.xdd(dt)
        for i in range(3):
            t = abs(acc1[i])
            if t > amax:
                amax = t
            t = abs(acc2[i])
            if t > amax:
                amax = t
        return amax


    #
    # it seems constrain_A can be identical for 3D and 1D(!)
    #
    def constrain_A(self):
        # hack for fixed dt
        #dt = DT_TEST
        #self.compute(dt)

        #
        # adaptive dt
        #
        dt = 0.2
        ni = 0
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.1
            else:
                break
        dt *= 0.9
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.02
            else:
                break
        print('Constrain_A: iterations: ', ni, ' dt = ', dt)
        print('             self.get_Amax(dt):',self.get_Amax(dt))
        self.dt = dt
        self.constrained = True
        self.compute(self.dt)
        return dt

    # also essentially unchanged!
    def timeEvolution(self):
        if not self.computed and not self.constrained:
            error('Cant compute timeEvolution until trajectory is computed and constrained')
        Np = 20-1
        x = []
        v = []
        a = []
        t = []
        for t1 in range (Np):
            t2 = self.dt*t1/Np
            t.append(t2)
            x.append(self.x(t2))
            v.append(self.xd(t2))
            a.append(self.xdd(t2))
        t.append(self.dt)
        x.append(self.x(self.dt))
        v.append(self.xd(self.dt))
        a.append(self.xdd(self.dt))
        ##print('A shape1: ',a.shape())
        t = np.array(t).T
        x = np.array(x).T
        v = np.array(v).T
        a = np.array(a).T
        print('X shape: (tev)',x.shape)

        return t,x,v,a

    #
    #   Energy cost: sum of squared acceleration
    def cost_e(self,a):
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is computed and constrained')
        c = 0.0
        n = len(a)
        for i in range(3):
            for a1 in a:
                c += self.dt*a1[i]*a1[i]/n
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        return self.dt


class trajectory1D:
    def __init__(self,p1,p2):
        print('trajectory1D: p1,p2: ',p1,p2)
        self.p1 = p1
        self.p2 = p2
        self.computed = False
        self.constrained = False
        self.valid = True

    def compute(self,dt):
        p1 = self.p1
        p2 = self.p2
        # p(0) = a0
        self.a0 = p1.x
        # v(0) = a1
        self.a1 = p1.v
        # p(dt) = a0 + self.a1*dt
        #define some constants
        dx = p2.x-p1.x
        dv = p2.v-p1.v
        b0 = dt
        b1 = dt*dt
        b2 = dt*dt*dt
        b3 = 2.0*dt
        b4 = 3.0*dt*dt
        # a3 solution:
        self.a3 = (dv - (b3/b1)*(dx-p1.v*b0))/(b4-b2*b3/b1)
        self.a2 = (dx - p1.v*b0 - self.a3*b2)/b1
        self.computed = True

    def x(self,t):
        return self.a0 + self.a1*t +    self.a2*t*t +     self.a3*t*t*t
    def xd(self,t):
        return           self.a1   + 2.0*self.a2*t  + 3.0*self.a3*t*t
    def xdd(self,t):
        return                       2.0*self.a2    + 6.0*self.a3*t

    def get_Amax(self,dt):
        if not self.computed:
            error('cant get_Amax() unless tr.computed is True')
        return max(abs(self.xdd(0)),abs(self.xdd(dt)))

    def constrain_A(self):
        # hack for fixed dt
        #dt = DT_TEST
        #self.compute(dt)

        #
        # adaptive dt
        #
        dt = 0.2
        ni = 0
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.1
            else:
                break
        dt *= 0.9
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.02
            else:
                break
        print('Constrain_A: iterations: ', ni, ' dt = ', dt)
        print('             self.get_Amax(dt):',self.get_Amax(dt))
        self.dt = dt
        self.constrained = True
        self.compute(self.dt)
        return dt

    def timeEvolution(self):
        if not self.computed and not self.constrained:
            error('Cant compute timeEvolution until trajectory is computed and constrained')
        Np = 20
        x = []
        v = []
        a = []
        t = []
        for t1 in range (Np):
            t2 = self.dt*t1/Np
            t.append(t2)
            x.append(self.x(t2))
            v.append(self.xd(t2))
            a.append(self.xdd(t2))
        t.append(self.dt)
        x.append(self.x(self.dt))
        v.append(self.xd(self.dt))
        a.append(self.xdd(self.dt))
        return t,x,v,a

    #
    #   Energy cost: sum of squared acceleration
    def cost_e(self,a):
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is computed and constrained')
        c = 0.0
        n = len(a)
        for a1 in a:
            c += self.dt*a1*a1/n
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        return self.dt

    def __repr__(self):
        return str(self.p1) + ' ---> ' + str(self.p2)

class Cm:
    def __init__(self):
        self.m = [[ 0 for x in range(M)] for y in range(M)]

    def fill(self, grid):
        print('starting fill...')
        for i1 in range(N):  # go through first points
            for j1 in range(N):
                for i2 in range(N):  # for all 2nd points
                    for j2 in range(N):
                        r = i1*N + j1
                        c = i2*N + j2
                        #print('r,c:',r,c)
                        #print('i1, j1, i2, j2:',i1,j1,i2,j2)
                        p1 = grid.gr[i1][j1]
                        p2 = grid.gr[i2][j2]
                        t = trajectory1D(p1,p2)
                        #print('computing cost: ',t)
                        if (p1.x == p2.x and p1.v == p2.v):
                            t.valid = False  # eliminate self transitions
                        else:
                            t.compute(DT_TEST)
                            t.constrain_A()
                        self.m[r][c] = t
        print('done with fill...')

    def __repr__(self):
        res = 'Cm cost matrix:\n'
        for i in range(M):
            for j in range(M):
                res += ' {:5.1f}'.format(self.m[i][j].cost_t)
            res += '\n'
        return res


class path:
    def __init__(self,grid,Cm):
        sr = startrow
        sc = startcol
        mark = [True for x in range(M)]
        count = 0
        mark[sr*N+sc] = False # mark our starting point
        self.Tcost = 0.0
        self.path = []
        crow = sr*N + sc  # starting point in cost matrix
        firstrow = crow
        while len(self.path) < M-1:
            print('looking for next path pt: row: ',crow)
            cmin = 99999999
            cminidx = 0
            found = False
            for ccol in range(M):
                if mark[ccol]:
                    if Cm.m[crow][ccol].valid:  # first attempt: pick first min cost traj.
                                                           # and elim self transitions
                        x,v,a,t = Cm.m[crow][ccol].timeEvolution()
                        if costtype == 'energy':
                            ccost = Cm.m[crow][ccol].cost_e(a)
                        elif costtype == 'time':
                            ccost = Cm.m[crow][ccol].cost_t()
                        else:
                            error('unknown cost type: '+costtype)
                        #print('  cost:',ccost)
                        if ccost < cmin:
                            #print('newmin: ', ccost, ccol)
                            cmin = ccost
                            cminidx = ccol
                            found = True
            if not found:
                error('I didnt find a minimum cost!!'+ str(count))
            if not mark[cminidx]:
                error('new path point is marked already')
            mark[cminidx] = False  # do not visit this point again
            if not Cm.m[crow][cminidx].valid:
                print('crow', crow, 'ccol: ', ccol)
                error('path: invalid new point')
            t = Cm.m[crow][cminidx]
            if crow != firstrow:
                pprev = self.path[-1].p2
                pcurr = t.p1
                if pprev != pcurr:
                    print('crow/firstrow: ', crow, firstrow)
                    print('adding traj: ', t)
                    error('path trajectories dont connect! '+str(crow))
            self.path.append(t)
            print('              adding traj to path: ', self.path[-1])
            crow = cminidx
            self.Tcost += cmin
            print('Total path cost ({:}) = {:8.2f}: '.format(costtype,self.Tcost))

    def check(self):
        if len(self.path) != M-1:
            error('wrong path length '+str(len(self.path)))
        i=0
        for t in self.path:
            if not t.valid:
                error('path contains invalid traj: '+str(i))
            if i > 0:
                if self.path[i-1].p2 != t.p1:
                    error('discontinuous path: '+str(i-1)+' '+str(i))
            i+=1

    def __repr__(self):
        r = '\n'
        i = 0
        for p in self.path:
            r += str(p)+' -> '
            i += 1
            if i%7 == 0:
                r += '\n'

        r += '\n total cost: ' + '{:7.1f}'.format(self.Tcost) + '\n'
        return r

    def compute_curves(self,idx):
        curvepts_x = []
        curvepts_v = []
        for i,tr in enumerate(self.path):
            if i == idx or idx < 0:
                if i == len(self.path):
                    break
                if not tr.valid:
                    error('compute_curves(): I shouldnt have found an invalid trajectory in path: '+str(i))
                if tr is None:
                    error('null traj: '+ str(i) + str(tr))
                if not tr.computed and not tr.constrained:
                    error('Cant plot until trajectory is computed and constrained '+ str(i) + str( p))
                dt = tr.dt
                for i in range(NPC):
                    t = dt*i/NPC
                    curvepts_x.append(tr.x(t))
                    curvepts_v.append(tr.xd(t))
                curvepts_x.append(tr.x(dt))
                curvepts_v.append(tr.xd(dt))
        return curvepts_x,curvepts_v


    def plot(self,idx):
        x_values = [point.p1.x for point in self.path]  # starting values
        y_values = [point.p1.v for point in self.path]
        x_values.append(self.path[-1].p2.x)
        y_values.append(self.path[-1].p2.v)

        print('plotting {:} xy values.'.format(len(x_values)))

        plt.plot(x_values, y_values, 'o')
        plt.xlabel('X')
        plt.ylabel('\dot{X}')
        plt.title('Path through Grid: minimize {:}  Amax = {:2.1f} '.format(costtype,AMAX))
        plt.grid(True)


        arrow_positions = np.array([x_values, y_values]).T
        arrow_directions = np.diff(arrow_positions, axis=0)

        # Plot arrows on the path at regular intervals
        arrow_interval = len(self.path) // (M-1) # Change 5 to adjust the arrow density
        arrow_positions = arrow_positions[:-1:arrow_interval]
        arrow_directions = arrow_directions[::arrow_interval]

        tr = self.path[0]
        startpt = plt.Circle((tr.p1.x,tr.p1.v),0.05,color='green')
        ax = plt.gca()
        ax.add_patch(startpt)

        plt.quiver(
            arrow_positions[:, 0],
            arrow_positions[:, 1],
            arrow_directions[:, 0],
            arrow_directions[:, 1],
            scale=1,
            scale_units='xy',
            angles='xy',
            width=0.005,
            color='red'
        )


        cx, cy = self.compute_curves(idx)
        plt.plot(cx,cy,color='blue')
        axlim = 2
        ax.set_xlim([-axlim,axlim])
        ax.set_ylim([-axlim,axlim])
        plt.show()



def tests():
    c1 = Cm()
    print (c1.m)

    print('Cm tests: ')
    print(c1)

    print('grid1D tests:')
    gt = grid1D(N)
    print(gt)

    print('cost tests')

    p1 = point1D(0,0)
    p2 = point1D(N,N)
    print ('{:5.1f}'.format(trajectory1D(p1,p2).cost_e()))
    p3 = point1D(1,1)
    print ('{:5.1f}'.format(trajectory1D(p1,p3).cost_e()))

    print('fill tests')
    c1.fill(gt)

    r = 2
    c = 2
    print('visualize cost from point ({:},{:})'.format(r,c))
    c1.visualize(gt,r,c)

    # try a path:
    p = path(gt,c1,r,c)
    print(p)


if __name__ ==  '__main__':
    print('main starting:')
    tests()

