#!/usr/bin/python3

import logging
import numpy as np
import matplotlib.pyplot as plt

def error(msg):
    print('Error: ')
    print(msg)
    quit()

N = 4
M = N*N  # number of distances

AMAX = 8  #  normalize for now.
N_CMP_PTS = 10 # points only for visualization

class grid2D:
    def __init__(self, N):
        self.N = N
        self.gr = []
        for i in range(self.N):
            row = []
            for j in range(self.N):
                row.append(point2D(i,j))
            self.gr.append(row)

    def __repr__(self):
        txt = ''
        for i in range(self.N):
            txt += '       '
            for j in range(self.N):
                txt +=  ' ({:5.2f},{:5.2f})'.format(self.gr[i][j].x, self.gr[i][j].v)
            txt += '\n'
        return txt

class point2D:
    def __init__(self,i,j):
        self.row = i
        self.col = j
        self.x =     2*j/(N-1) - 1
        self.v = -1*(2*i/(N-1) - 1)
        self.tr = None

    def __repr__(self):
        return ' ({:4.1f}, {:4.1f})'.format(self.x, self.v)
        #return ' ({:}, {:})'.format(self.row, self.col)


class trajectory2D:
    def __init__(self,p1,p2):
        print('trajectory2D: p1,p2: ',p1,p2)
        self.p1 = p1
        self.p2 = p2
        self.computed = False
        self.constrained = False

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
        return max(abs(self.xdd(0)),abs(self.xdd(dt)))

    def constrain_A(self):
        dt = 2.0
        #dt = 0.2
        #ni = 0
        #while True:
            #ni += 1
            #self.compute(dt)
            #if self.get_Amax(dt) > AMAX:
                #dt *= 1.1
            #else:
                #break
        #print('Constrain_A: iterations: ', ni, ' dt = ', dt)
        #print('             self.get_Amax(dt):',self.get_Amax(dt))
        self.dt = dt
        self.constrained = True
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
    def cost_e(self,a,dt):
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is computed and constrained')
        c = 0.0
        n = len(a)
        for a1 in a:
            c += dt*a1*a1/n
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        return self.dt


    #def compute(self):
        #tij_x = []
        #tij_v = []
        #x = self.p1.x
        #v = self.p1.v
        #dx1 = self.dx/N_CMP_PTS
        #dv1 = self.dv/N_CMP_PTS
        #for i in range(N_CMP_PTS):
            #tij_x += dx1
            #tij_v += dv1
            #tij_x.append(x)
            #tij_v.append(v)
        #return( [tij_x, tij_v] )

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
                        t = trajectory2D(grid.gr[i1][j1], grid.gr[i2][j2])
                        t.compute(1.0)
                        t.constrain_A()
                        self.m[r][c] = t
        print('done with fill...')

    def visualize(self,grid,i1,j1):  # must be filled first
        p1 = grid.gr[i1][j1]
        index = i1*N+j1
        g = [[0 for x in range(N)] for y in range(N)]
        for i in range(N):
            for j in range(N):
                g[i][j] = self.m[index][i*N+j]
        # print it
        for i in range(N):
            for j in range(N):
                print(' {:4.1f}'.format(g[i][j]),end='')
            print('')

    def __repr__(self):
        res = 'Cm cost matrix:\n'
        for i in range(M):
            for j in range(M):
                res += ' {:5.1f}'.format(self.m[i][j])
            res += '\n'
        return res


class path:
    def __init__(self,grid,Cm,sr,sc):
        mark = [True for x in range(M)]
        mark[sr*N+sc] = False # mark our starting point
        self.Tcost = 0.0
        self.path = [ grid.gr[sr][sc] ]
        crow = sr*N + sc  # starting point in cost matrix
        while len(self.path) < M:
            cmin = 999999
            cminidx = 0
            print('looking for next path pt: ')
            for ccol in range(M):
                if mark[ccol]:  # first attempt: pick first min cost traj.
                    Cm.m[crow][ccol].compute(1.0)
                    Cm.m[crow][ccol].constrain_A()
                    #ccost = Cm.m[crow][ccol].cost_e
                    ccost = Cm.m[crow][ccol].cost_t()
                    print('  cost:',ccost)
                    if ccost < cmin and ccost > 0.0:
                        print('newmin: ', ccost)
                        cmin = ccost
                        cminidx = ccol
            mark[cminidx] = False
            mingrow = cminidx//N
            mingcol = cminidx-mingrow*N
            self.path.append(grid.gr[mingrow][mingcol])
            self.path[-1].tr = Cm.m[crow][ccol]
            print('              adding pt: ', self.path[-1])
            crow = cminidx
            self.Tcost += cmin
            cminidx += 1
        # now put in trajectory for starting point
        self.path[0].tr = trajectory2D(grid.gr[sr][sc], self.path[1])
        self.path[0].tr.compute(1.0)
        self.path[0].tr.constrain_A()

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

    def cc2(self,idx):
        curvepts_x = []
        curvepts_v = []
        npc = 20
        for i,p in enumerate(self.path):
            if i == idx:
                if i+1 == len(self.path):
                    break
                tr = p.tr
                if tr is None:
                    error('null traj: '+ str(i) + str( p))
                if not tr.computed and not tr.constrained:
                    error('Cant plot until trajectory is computed and constrained '+ str(i) + str( p))
                dt = tr.dt
                for i in range(npc):
                    t = dt*i/npc
                    curvepts_x.append(tr.x(t))
                    curvepts_v.append(tr.xd(t))
                curvepts_x.append(tr.x(dt))
                curvepts_v.append(tr.xd(dt))
        return curvepts_x,curvepts_v

    #def compute_curves(self):
        #curvepts_x = []
        #curvepts_v = []
        #npc = 12
        #for i,p in enumerate(self.path):
            #if i+1 == len(self.path):
                #break
            #dx = self.path[i+1].x - p.x
            #dv = self.path[i+1].v - p.v
            #dt = abs(dv/AMAX)
            #if dt > 0:
                #a = 0.5*dv/dt
            #else:
                #a = 0.0
                #dt = dx/p.v
            #for i in range(npc):
                #t = dt*i/npc
                #curvepts_x.append(p.x + p.v*t + 0.5*a*t*t)
                #curvepts_v.append(      p.v   +     a*t)
        #return curvepts_x,curvepts_v

    def plot(self,idx):
        x_values = [point.x for point in self.path]
        y_values = [point.v for point in self.path]

        print('plotting {:} xy values.'.format(len(x_values)))

        plt.plot(x_values, y_values, '-o')
        plt.xlabel('X')
        plt.ylabel('\dot{X}')
        plt.title('Path through Grid')
        plt.grid(True)


        arrow_positions = np.array([x_values, y_values]).T
        arrow_directions = np.diff(arrow_positions, axis=0)

        # Plot arrows on the path at regular intervals
        arrow_interval = len(self.path) // 16 # Change 5 to adjust the arrow density
        arrow_positions = arrow_positions[:-1:arrow_interval]
        arrow_directions = arrow_directions[::arrow_interval]

        pt = self.path[0]
        startpt = plt.Circle((pt.x,pt.v),0.05,color='green')
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


        cx, cy = self.cc2(idx)
        plt.plot(cx,cy,color='blue')
        ax.set_xlim([-3,3])
        ax.set_ylim([-3,3])
        plt.show()



def tests():
    c1 = Cm()
    print (c1.m)

    print('Cm tests: ')
    print(c1)

    print('grid2D tests:')
    gt = grid2D(N)
    print(gt)

    print('cost tests')

    p1 = point2D(0,0)
    p2 = point2D(N,N)
    print ('{:5.1f}'.format(trajectory2D(p1,p2).cost_e()))
    p3 = point2D(1,1)
    print ('{:5.1f}'.format(trajectory2D(p1,p3).cost_e()))

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

