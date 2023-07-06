#!/usr/bin/python3

import numpy as np
import math
import matplotlib.pyplot as plt
import itertools as itt
import datetime
import brl_data.brl_data as bd
import random

def error(msg):
    print('Error: ')
    print(msg)
    quit()

N = 4

AMAX = 2  #  normalize for now.
DT_TEST = 2.0


NPC = 20  # number plot points per curve

costtype = 'time'
costtype = 'energy'

startrow = 3
startcol = 1

def configure(fp=None):
    global costtype, AMAX, DT_TEST, N, NPC, startrow, startcol
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

def ij2idx(i,j):
    return i*N+j
def idx2ij(idx):
    i = idx//N
    j = idx-N*i
    return i,j


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
        if i > N or j>N or i<=-1 or j<=-1:
            error('point2D i,j is too big: '+str(i)+' '+str(j)+ ' /'+str(N))
        self.row = i
        self.col = j
        self.x =     2*j/(N-1) - 1
        self.v = -1*(2*i/(N-1) - 1)
        if self.x < -1.05 or self.x > 1.05:
            msg = 'point2D: Im creating bogus point coordinates: '
            msg += str(i) + '  '
            msg += str(j) + '  '
            msg += str(self.x) + '  '
            msg += str(self.v) + '  '
            error(msg)
        self.tr = None

    def __eq__(self,x):
        return self.row == x.row and self.col == x.col

    def __repr__(self):
        return ' ({:4.1f}, {:4.1f})'.format(self.x, self.v)
        #return ' ({:}, {:})'.format(self.row, self.col)


class trajectory2D:
    def __init__(self,p1,p2):
        #print('trajectory2D: p1,p2: ',p1,p2)
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
        #print('Constrain_A: iterations: ', ni, ' dt = ', dt)
        #print('             self.get_Amax(dt):',self.get_Amax(dt))
        self.dt = dt
        self.constrained = True
        self.compute(self.dt) # because we changed dt
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

class Cm: # matrix full of trajectory2D objs
    def __init__(self):
        self.m = [[ 0 for x in range(N*N)] for y in range(N*N)]

    def fill(self, grid):
        print('starting fill...{:}x{:}'.format(N,N))
        for i1 in range(N):  # go through first points
            for j1 in range(N):
                for i2 in range(N):  # for all 2nd points
                    for j2 in range(N):
                        m=max(i1,j1,i2,j2)
                        if m>N-1:
                            error('index too big: {:}/{:}'.format(m,N))
                        r = i1*N + j1
                        c = i2*N + j2
                        #print('r,c:',r,c)
                        #print('i1, j1, i2, j2:',i1,j1,i2,j2)
                        p1 = grid.gr[i1][j1]
                        p2 = grid.gr[i2][j2]
                        t = trajectory2D(p1,p2)
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
        for i in range(N*N):
            for j in range(N*N):
                res += ' {:5.1f}'.format(self.m[i][j].cost_t)
            res += '\n'
        return res

class quartile():
    def __init__(self):
        self.qset = set()
        self.n = 0

    def add(self,c):
        self.qset.add(c)
        self.n += 1

    def purge(self,NptsQtile,sbelow, sabove):
        change = False
        big = False
        if self.n > NptsQtile:
            big = True
            smax = -1
            smin = 9999999999999
            for se in self.qset:
                if  se > smax:
                    smax = se
                if  se < smin:
                    smin = se
        if big and sbelow:
            sbelow.qset.add(smin)
            self.qset.remove(smin)
            self.n -= 1
            change = True
        if big and sabove:
            if smax == -1:
                error('somethings wrong with quartile '+str(self.n)+', '+str(self.qset))
            sabove.qset.add(smax)
            self.qset.remove(smax)
            self.n -= 1
            change = True
        return change

class path:
    def __init__(self,grid,Cm):
        self.Cm = Cm      # cost matrix (actually trajectories)
        self.grid = grid
        self.sr = startrow
        self.sc = startcol
        self.mark = [True for x in range(N*N)]
        count = 0
        self.mark[self.sr*N+self.sc] = False # mark our starting point
        self.Tcost = 0.0
        self.path = []

    def sampleSearch(self,Cm,ns=1000):
        return self.bruteForce(Cm,sampling=True,n=ns)

    def bruteForce(self,Cm,sampling=False,n=0): # path class
        n_all_paths = math.factorial(N*N)
        print('Starting brute force: N=',N)
        if sampling:
            print('   Sampling {:} paths out of {:12.5e}'.format(n,float(n_all_paths)))

        LOWMEM = False
        if N>3:
            LOWMEM = True

        SPEEDTEST = False  # just run 2000 paths to measure speed

        STOREDATA = True   # write all perms and costs to a data file

        if STOREDATA:
            dfbf = bd.datafile('TSP_perms','BH','simulation')
            dataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
            codeFolder = ''
            dfbf.set_folders(dataFolder,codeFolder)
            print('Saving permutations (paths) to: ',dfbf.name)
            itype = str(type(5))
            ftype = str(type(3.1415))
            tps = [itype]*(N*N)      # path point seq
            names = [' ']*(N*N)
            tps.append(ftype) # the path cost's type
            names.append('Cost')
            dfbf.metadata.d['Types'] = tps
            dfbf.metadata.d['Names'] = names
            dfbf.metadata.d['CostType'] = costtype
            #
            dfbf.open()  # let's open the file (default is for writing)

        ##1) list all possible paths
        if not sampling:
            print('We are about to find all paths through ',N*N,' nodes')
            x=input('ready?..')
            piter = itt.permutations(range(N*N),N*N) # not a list!
        else:
            print('We are generating {:} random paths through {:} nodes'.format(n,N*N))
            m = N*N
            piter = []
            phashlist = []
            for i in range(n):
                p = list(range(N*N))
                random.shuffle(p)
                pthash = 0
                for j in p:
                    pthash += j
                    pthash *= N*N
                if pthash not in phashlist: # we've found a new pt
                    piter.append(p)
                    phashlist.append(pthash)

        print('Path enumeration complete:')
        secPerLoop = 0.0003366 # measured on IntelNUC
        secPerLoop = 0.0008419 # Dell XPS-13
        print('   Estimated completion time: ',n_all_paths*secPerLoop,' sec.')
        hrs = n_all_paths*secPerLoop/(60*60)
        print('   Estimated completion time: ',hrs,' hrs.')
        days = hrs/24
        print('   Estimated completion time: ',days,' days.')
        months = 12*days/365
        print('   Estimated completion time: ',months,' months.')
        years = months/12
        print('   Estimated completion time: ',years,' years.')

        x=input('ready?..')

        #2) evaluate their costs
        path_costs = []
        n = -1
        cmin = 99999999999
        pmin = []
        if SPEEDTEST:
            Navg = 20000
            ts1 = datetime.datetime.now()
        for p in piter:  # piter returns set of points
            p = list(p)
            n += 1
            #print('path: ',n, p)
            c = 0.0
            tmpPath = []
            for i in range(len(p)-1):
                # build next trajectory
                row = p[i]//N
                col = p[i]-row*N
                p1 = point2D(row,col)
                row = p[i+1]//N
                col = p[i+1]-row*N
                p2 = point2D(row,col)
                tr = trajectory2D(p1,p2)
                tr.constrain_A()
                tmpPath.append(tr)
                if costtype == 'energy':
                    x,v,a,t = tr.timeEvolution()
                    c += tr.cost_e(a)
                elif costtype == 'time':
                    c += tr.cost_t()
                else:
                    error('unknown cost type: '+costtype)
            if SPEEDTEST and n > Navg:
                break
            if n%2000 == 0:
                    print('path ',n)
            if STOREDATA:
                row = p
                row.append(c)
                dfbf.write(row)
            if c < cmin:
                cmin = c
                pmin = path(self.grid,self.Cm)
                pmin.path = tmpPath
                pmin.Tcost = c
                nmin = n
            #print(' path cost: {:4.2f}'.format(c))
            if not LOWMEM:
                path_costs.append(c)
        if SPEEDTEST:
            print('timed {:} paths.'.format(Navg))
            ts2 = datetime.datetime.now()
            dt = (ts2-ts1).total_seconds()
            print('seconds per {:} paths: {:}'.format(Navg, float(dt)))
            print('seconds per path: {:}'.format(float(dt)/Navg))
            quit()
        if not LOWMEM:
            print('quantiling the costs')
            q = [0.25,0.5,0.75,1.0]
            qs = np.quantile(path_costs, q )
            minCost= np.min(path_costs)
            maxCost= np.max(path_costs)
            print('Quartile Report:')
            print('  Min:        {:4.1f}'.format( minCost))
            print('    Q1: {:4.2f}  {:4.1f}'.format( q[0], qs[0]))
            print('    Q2: {:4.2f}  {:4.1f}'.format( q[1], qs[1]))
            print('    Q3: {:4.2f}  {:4.1f}'.format( q[2], qs[2]))
            print('    Q4: {:4.2f}  {:4.1f}'.format( q[3], qs[3]))
            print('  Max:        {:4.1f}'.format( maxCost ))
        else: # not enough memory for quartiles(!)
            print('Lowest cost path: ', p)
            print(pmin, 'path cost: ', cmin)

        if STOREDATA:
            dfbf.metadata.d['Min Cost']=minCost
            dfbf.metadata.d['Max Cost']=maxCost
            dfbf.metadata.d['Quartiles']=list(qs)
            dfbf.close()
        return pmin, pmin.Tcost

    def heuristicSearch(self):  # path class
        # add to self.path[] one traj at a time greedily
        crow = self.sr*N + self.sc  # starting point in cost matrix
        firstrow = self.sr*N + self.sc
        while len(self.path) < N*N-1: # build path up one pt at a time
            print('looking for next path pt: row: ',crow)
            cmin = 99999999
            cminidx = 0
            found = False
            for ccol in range(N*N):
                if self.mark[ccol]:
                    if self.Cm.m[crow][ccol].valid:  # first attempt: pick first min cost traj.
                                                           # and elim self transitions
                        x,v,a,t = self.Cm.m[crow][ccol].timeEvolution()
                        if costtype == 'energy':
                            ccost = self.Cm.m[crow][ccol].cost_e(a)
                        elif costtype == 'time':
                            ccost = self.Cm.m[crow][ccol].cost_t()
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
            if not self.mark[cminidx]:
                error('new path point is marked already')
            self.mark[cminidx] = False  # do not visit this point again
            if not self.Cm.m[crow][cminidx].valid:
                print('crow', crow, 'ccol: ', ccol)
                error('path: invalid new point')
            t = self.Cm.m[crow][cminidx]  # get current best trajectory
            if crow != firstrow:  # no trajectory has t.p2=startPoint
                pprev = self.path[-1].p2
                pcurr = t.p1
                if pprev != pcurr: # check for error
                    print('crow/firstrow: ', crow, firstrow)
                    print('adding traj: ', t)
                    error('path trajectories dont connect! '+str(crow))
            self.path.append(t)
            print('              adding traj to path: ', self.path[-1])
            crow = cminidx
            self.Tcost += cmin
            print('Total path cost ({:}) = {:8.2f}: '.format(costtype,self.Tcost))
        return self, self.Tcost

    def check(self):
        if len(self.path) != N*N-1:
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
        if type(idx) == type(5):
            idx = [idx]
        elif type(idx) != type([5]):
            error('compute_curves: idx is not an int or a list')
        curvepts_x = []
        curvepts_v = []
        if idx[0] < 0: # idx=-1 is a flag for all of the path
            idx = range(len(self.path))
        for i,tr in enumerate(self.path):
            if i in idx:
                if i == len(self.path):
                    break
                if not tr.valid:
                    error('compute_curves(): I shouldnt have found an invalid trajectory in path: '+str(i))
                if tr is None:
                    error('null traj: '+ str(i) + str(tr))
                if not tr.computed and not tr.constrained:
                    error('Cant plot until trajectory is computed and constrained '+ str(i) + str(tr))
                dt = tr.dt
                for i in range(NPC):
                    t = dt*i/NPC
                    curvepts_x.append(tr.x(t))
                    curvepts_v.append(tr.xd(t))
                curvepts_x.append(tr.x(dt))
                curvepts_v.append(tr.xd(dt))
        return curvepts_x,curvepts_v


    def plot(self,idx, note=''): # plot a path with trajectories
        self.plotOnePath(self.plotSetup(note))
        self.plotDone()

    def plotSetup(self,note):
        fig = plt.figure()
        plt.title('Path through Grid: minimize {:}  Amax = {:2.1f}\n            {:} '.format(costtype,AMAX,note))
        plt.xlabel('X')
        plt.ylabel('\dot{X}')
        plt.grid(True)
        return fig

    def plotOnePath(self,fig):
        x_values = [point.p1.x for point in self.path]  # starting values
        y_values = [point.p1.v for point in self.path]
        x_values.append(self.path[-1].p2.x)
        y_values.append(self.path[-1].p2.v)
        ax = plt.gca()
        ax.plot(x_values, y_values, 'o')
        arrow_positions = np.array([x_values, y_values]).T
        arrow_directions = np.diff(arrow_positions, axis=0)

        # Plot arrows on the path at regular intervals
        arrow_interval = len(self.path) // (N*N-1) # Change 5 to adjust the arrow density
        #print('arrow_interval: {:}  len(self.path) {:}  N*N {:}'.format(arrow_interval, len(self.path), N*N))
        arrow_positions = arrow_positions[:-1:arrow_interval]
        arrow_directions = arrow_directions[::arrow_interval]

        tr = self.path[0]
        startpt = plt.Circle((tr.p1.x,tr.p1.v),0.05,color='green')
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


        cx, cy = self.compute_curves(-1) #compute trajectory path
        ax.plot(cx,cy,color='blue')
        axlim = 2
        ax.set_xlim([-axlim,axlim])
        ax.set_ylim([-axlim,axlim])

    def plotDone(self):
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
    p.heuristicSearch()

    print(p)


if __name__ ==  '__main__':
    print('main starting:')
    tests()

