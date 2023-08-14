#!/usr/bin/python3

import os
import json
import numpy as np
import math
import matplotlib.pyplot as plt
import itertools as itt
import datetime
import brl_data.brl_data as bd
import random
import socket
import time


def error(msg):
    print('Error: ')
    print(msg)
    quit()

PCNAME = str(socket.gethostname())
N = 4

AMAX = 2  #  normalize for now.
DT_TEST = 2.0


NPC = 20  # number plot points per curve

costtype = 'time'
costtype = 'energy'

startrow = 3
startcol = 1

def configure(fp=None):
    global costtype, AMAX, DT_TEST, N, M, NPC, startrow, startcol
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

#  'i,j' 'r,c', 'row,col' are used kind of interchangeably
#   to indicate point in rectangular grid ... sorry!
def ij2idx(i,j):   #  0-N*N-1
    return i*N+j
def idx2ij(idx):   #  0-N, 0-N
    i = idx//N
    j = idx-N*i
    return i,j

def pt2idx(pt):
    return ij2idx(pt.row,pt.col)

def idx2rc(idx):
    # index to point row, col
    # same as above
    return idx2ij(idx)

def predict_timing(df, searchtype, nsamp):
    #
    # predict the search timing
    #
    savedir = ''
    if df is not None:
        savedir = df.folder
    filename = savedir + 'searchTiming.json'
    OK = True
    key = searchtype + '-' + PCNAME
    if os.path.isfile(filename):
        fd = open(filename,'r')
        d = json.load(fd)
        try:
            t = d[key]
        except:
            OK=False
        if OK:
            print('your predicted search time is:')
            print('search type/PC:',key)
            print('rate:          ',d[key],'/sec')
            sec = float(d[key])*nsamp
            mins = sec/60
            hrs = mins/60
            days = hrs/24
            years = days/365
            print('predicted time: ')
            fmth='{:30} {:>12} {:>12} {:>12} {:>12}'
            print(fmth.format('PC type','mins','hrs','days','years'))
            fmts='{:30} {:12.2f} {:12.2f} {:12.2f} {:12.2f}'
            print(fmts.format(key, mins,hrs,days,years))
            if mins>2:
                x=input('OK to continue? ...')
    else:
        OK=False #path is not a file
    if not OK:
        print('no search speed info available for your configuration: ',key)
        x=input('OK to continue? ...')

def save_timing(df, searchname,systemName,rate):
    savedir = ''
    if df is not None:
        savedir = df.folder
    filename = savedir + 'searchTiming.json'
    if os.path.isfile(filename):
        fd = open(filename,'r')
        d = json.load(fd)
    else:
        d = {}
    # add to the dict and resave
    d[searchname+'-'+systemName] = rate
    fd = open(filename,'w')
    json.dump(d,fd,indent=4)
    return


def plotSave(fig, dpi, imagedir, imagename):
    #fig = plt.gcf()
    #datad = string path to data directory
    if imagedir[-1] != '/':
        imagedir += '/'
    imagepath = imagedir+imagename
    if '.png' not in imagepath:
        imagepath += '.png'
    fig.savefig(imagepath,dpi=dpi)
    print('your plot is saved to: ',imagepath)

    ####  keep a "log book" if saved
    dim = '6D'
    notes = '{:}, {:}'.format(dim,imagepath)
    now = datetime.datetime.now()
    dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
    logentry = '{:}, plot saved: {:}'.format(dtstring,notes)
    fdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    f =open(fdir + 'image_log.txt','a')
    print(logentry,file=f)
    f.close()
    ####


class grid2D:
    def __init__(self, N):
        self.N = N
        self.gr = []
        self.randgrid = False
        self.fromFile = False
        for i in range(self.N):
            row = []
            for j in range(self.N):
                row.append(point2D(i,j))
            self.gr.append(row)

    def savePoints2D(self,df):   # save points generated into a file
        if not self.randgrid:
            error('savePoints2D: Do not save points if they are not random')
        # each entry in Cm is a trajectory2D  we need only store p1 and p2
        it =str(type(5))
        ft =str(type(3.14159))
        df.metadata.d['Types'] = [it,it,ft,ft]
        df.metadata.d['Names'] = ['i1','j1','x1','v1']
        df.metadata.d['Ncols'] = len(df.metadata.d['Names'])
        df.descript_str = 'randomGridPointSet' # standardize for point storage filename fields
        df.metadata.d['Research Question'] = 'RandomGridPointSet' # standardize for RQ
        df.open('w')
        # write out key info for each point.
        for i in range(N):
            for j in range(N):
                p1 = self.gr[i][j]
                row = [i,j, p1.x, p1.v ]
                df.write(row)
        df.close()

    def readPoints2D(self,df):  # read randomized points from a file
        if not self.randgrid:
            error('readPoints2D: Do not read in points if they are not random')
        df.open('r')
        # read in pair of points for each Cm.m entry and store trajectory.
        ptindex = 0
        for row in df.reader:
            i = int(row[0])  # 0 -- N*N-1
            j = int(row[1])
            print(f"I'm creating point at {i}, {j}")
            p1 = point2D(i,j)
            p1.x = float(row[2])
            p1.v = float(row[3])
            self.gr[i][j]=p1
            ptindex += 1
        df.close()
        self.fromFile = True
        # return the file name from which the pts were read for the record
        return df.hashcode

    def __repr__(self):
        txt = ''
        for i in range(self.N):
            txt += '       '
            for j in range(self.N):
                txt +=  ' ({:5.2f},{:5.2f})'.format(self.gr[i][j].x, self.gr[i][j].v)
            txt += '\n'
        return txt

class point2D:
    def __init__(self,i,j):  # index the spatial grid 0--N-1, 0--N-1
        if i > N or j>N or i<=-1 or j<=-1:
            error('point2D i,j is too big: '+str(i)+' '+str(j)+ ' /'+str(N))
        self.row = i
        self.col = j
        self.x =     2*j/(N-1) - 1
        self.v = -1*(2*i/(N-1) - 1)
        if self.x < -1.05 or self.x > 1.05:
            msg = f"point2D: I'm creating bogus point coordinates: {i} {j} {self.x} {self.v}"
            error(msg)
        self.tr = None

    def randomize(self):
        self.x = np.random.uniform(-1,1)
        self.v = np.random.uniform(-1,1)
        return

    def __eq__(self,x):
        #return self.row == x.row and self.col == x.col
        epsilon = 0.002
        ex = abs(self.x-x.x)/x.x < epsilon
        ev = abs(self.v-x.v)/x.v < epsilon
        return ex and ev

    def __repr__(self):
        return ' ({:4.2f}, {:4.2f})'.format(self.x, self.v)
        #return ' ({:}, {:})'.format(self.row, self.col)


class trajectory2D:
    def __init__(self,p1,p2):
        #print('trajectory2D: p1,p2: ',p1,p2)
        self.p1 = p1
        self.p2 = p2
        self.e_cost = None
        self.t_cost = None
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

    def timeEvolution(self,ACC_ONLY=False):
        if not self.computed and not self.constrained:
            error('Cant compute timeEvolution until trajectory is computed and constrained')
        Np = 20
        x = []
        v = []
        a = []
        t = []
        if not ACC_ONLY:
            for t1 in range (Np):
                t2 = self.dt*t1/Np
                t.append(t2)
                x.append(self.x(t2))
                v.append(self.xd(t2))
                a.append(self.xdd(t2))
            # add last points to close
            t.append(self.dt)
            x.append(self.x(self.dt))
            v.append(self.xd(self.dt))
            a.append(self.xdd(self.dt))
            return t,x,v,a
        if ACC_ONLY:
            for t1 in range (Np):
                t2 = self.dt*t1/Np
                a.append(self.xdd(t2))
            a.append(self.xdd(self.dt))
            return a

    #   Energy cost: sum of squared acceleration
    def cost_e(self,a):
        #
        # Note: there is an analytical expression for this:
        #  \int_0^{\delta t} \ddot{x} = 4a_2(\delta t) + 12a_2a_3(\delta t)^2 + 12a_3^2(\delta t)^3
        #    (derived 8/14/23)
        #
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is computed and constrained')
        c = 0.0
        n = len(a)
        for a1 in a:
            c += self.dt*a1*a1/n
        self.e_cost = c
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        self.t_cost = self.dt
        return self.dt

    def __repr__(self):
        return str(self.p1) + ' ---> ' + str(self.p2)

    def aprtr(self,msg):
        print(msg)
        print('start: {:4.2f}, {:4.2f}'.format(self.p1.x,self.p1.v))
        print('goal:  {:4.2f}, {:4.2f}'.format(self.p2.x,self.p2.v))

class Cm: # matrix full of trajectory2D objs
    def __init__(self,df=None):
        self.m = [[ 0 for x in range(N*N)] for y in range(N*N)]
        # normally the points are in a regular 6 dim grid.  if Randgrid
        # is true, they will be converted into uniform([-1,1)) in all coordinates
        self.randgrid = False
        if df is not None:
            df.metadata.d['Random Grid'] = False


    def set_GridRandomize(self,df=None):
        self.randgrid = True
        if df is not None:
            df.metadata.d['Random Grid'] = True
        return

    def fill(self, grid):
        nselftr = 0
        print('starting fill...{:}x{:}'.format(N,N))
        grid.randgrid = False # we do the setting of this and the randomization itself HERE
        if grid.fromFile:  # we have read grid already from file and it must be random
            grid.randgrid=True
        elif self.randgrid:  # but we haven't read a file (gen from scratch)
            grid.randgrid=True
            # go through the grid and randomize their x,v
            for r in range(N):  # for all 2nd points
                for c in range(N):
                    p1 = grid.gr[r][c]
                    p1.randomize()
                    grid.gr[r][c] = p1

        # now get a tr (and costs) for all possible transitions (row x col)
        for i1 in range(N*N):
            for j1 in range(N*N):
                r1,c1 = idx2ij(i1)
                r2,c2 = idx2ij(j1)
                p1 = grid.gr[r1][c1]
                p2 = grid.gr[r2][c2]
                # we now have p1 and p2
                t = trajectory2D(p1,p2)
                if (p1.x == p2.x and p1.v == p2.v): # robust to random pts
                    nselftr += 1
                    t.valid = False  # eliminate self transitions
                else:  # compute the two costsF
                    #t.compute(DT_TEST) #get coeffs
                    t.constrain_A()    #constrain for Amax
                    a = t.timeEvolution(ACC_ONLY=True)
                    #compute traj costs
                    t.cost_e(a)
                    t.cost_t()
                self.m[i1][j1] = t   # store the trajectory
        print('done with fill...')
        print('# of self-state (invalid) transitions: ',nselftr, ' (expected N*N): ', N*N)

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

MAXTIEHISTO = 80  # how many bins for the tie frequency histogram

class path:
    def __init__(self,grid,Cm):
        self.Cm = Cm      # cost matrix (actually trajectories)
        self.grid = grid
        self.sr = startrow
        self.sc = startcol
        self.mark = [True for x in range(N*N)]
        self.mark[self.sr*N+self.sc] = False # mark our starting point
        self.Tcost = 0.0
        self.path = []  # the path as a list of trajectories
        self.idxpath = [] # the path as a list of indeces (0..N*N)
        self.searchtype = 'none yet'
        self.datafile = None
        #collect tie stats on this path
        self.maxTiesHSearch = -99999999  # most ties when greedy searching
        self.tie_freq = np.zeros(MAXTIEHISTO)  # histogram of how many ties of
        self.costset = set()  # testing only

    def search(self,searchtype,dfile=None,nsamples=1000):
        self.searchtype = searchtype

        predict_timing(dfile, searchtype, nsamples)
        #
        #  start timer
        ts1 = datetime.datetime.now()
        #
        #  select the type of search to do
        #
        if searchtype.startswith('heur'):
            p, cmin = self.heuristicSearch()
        elif searchtype.startswith('multi'):
            if dfile is None:
                error('path.search: multi heuristic search requires a dfile')
            p, cmin = self.multiHSearch(dfile,nsamples)
        elif searchtype.startswith('exh'):
            if dfile is None:
                error('path.search: brute force (exhaustive) search requires a dfile')
            p, cmin = self.bruteForce(dfile=dfile)
        elif searchtype.startswith('sampling'):
            if dfile is None:
                error('path.search: sampling search requires a dfile')
            p, cmin = self.sampleSearch(dfile=dfile,nsamples=nsamples)
        else:
            error('path.search: unknown search type: ', searchtype)
        #report timing
        ts2 = datetime.datetime.now()
        dt = (ts2-ts1).total_seconds()
        print('seconds per {:} paths: {:}'.format(nsamples, float(dt)))
        print('seconds per path: {:}'.format(float(dt)/nsamples))
        save_timing(dfile,searchtype,PCNAME,float(dt)/nsamples)
        p.datafile=dfile
        return p,cmin

    def cost(self):  # path.cost()
        c = 0.0
        if costtype == 'time':
            for t in self.path:
                c += t.t_cost
        elif costtype == 'energy':
            for t in self.path:
                c += t.e_cost
        else:
            error('path.cost() illegal cost type: ',costtype)
        print('                                       tcost: ',c)
        return c

    def sampleSearch(self,dfile=None,nsamples=977):
        return self.bruteForce(dfile=dfile,sampling=True,nsamples=nsamples)

    def bruteForce(self,dfile=None,sampling=False,nsamples=0): # path class
        if self.searchtype.startswith('none'):
            self.searchtype = 'exhaustive'
        n_all_paths = math.factorial(N*N)
        print('Starting {:} search: N={:}'.format(self.searchtype,N))
        if sampling:
            print('   Sampling {:} paths out of {:12.5e}'.format(nsamples,float(n_all_paths)))

        LOWMEM = False
        if N>3:
            LOWMEM = True  # we just can't store that many paths

        if dfile is not None:
            self.datafile = dfile
            STOREDATA = True   # write all perms and costs to a data file
        else:
            STOREDATA = False

        if STOREDATA:
            dfbf = dfile
            print('Saving permutations (paths) to: ',dfbf.name)
            itype = str(type(5))
            ftype = str(type(3.1415))
            tps = [itype]*(N*N)      # path point seq
            tps.append(ftype) # the path cost's type
            names = []
            for i in range(N*N):
                names.append(f'p{i}')
            names.append('Cost')
            dfbf.metadata.d['Ncols'] = len(names)
            dfbf.metadata.d['Types'] = tps
            dfbf.metadata.d['Names'] = names
            dfbf.metadata.d['CostType'] = costtype
            dfbf.metadata.d['SearchType'] = self.searchtype
            #dfbf.metadata.d['#samples'] = nsamples  SHOULD BE AUTO

            #
            dfbf.open()  # let's open the file (default is for writing)

        ##1) list all possible paths
        if not sampling:
            print('We are about to find all paths through ',N*N,' nodes')
            x=input('ready?..')
            piter = itt.permutations(range(N*N),N*N) # not a list!
        else:
            print('We are generating {:} random paths through {:} nodes'.format(nsamples,N*N))
            piter = []
            phashset = set()
            for i in range(nsamples):
                p = list(range(N*N))
                random.shuffle(p) # generate a path as random list of indices
                pthash = ''
                for j in p:
                    pthash +='{:3d}'.format(j) # 'hash' the list
                if pthash not in phashset: # we've found a new pt
                    piter.append(p)
                    phashset.add(pthash)

        print('Path enumeration complete:')

        # keep around for history
        #secPerLoop = 0.0003366 # measured on IntelNUC
        #secPerLoop = 0.0008419 # Dell XPS-13

        #2) evaluate their costs
        path_costs = []
        n = -1
        cmin = 99999999999
        cmax = 0
        pmax = path(self.grid,self.Cm)   # storage for int results
        pmin = path(self.grid, self.Cm)

        itrct = -1
        for p in piter:  # piter returns a whole "list" of point indices each time
            itrct += 1
            idxpath = list(p)
            p = path(self.grid, self.Cm) # temp variable
            p.idxpath = idxpath # use own idxpath for tmp idx path storage(!)
            p.path = []  # use own path for tmp path storage(!)
            iterct = -1
            for i,ptidx in enumerate(p.idxpath[1:]): # list of points
                iterct+=1
                prvidx = p.idxpath[i]  # kind of actually "i-1" - b/c [1:]
                tr = self.Cm.m[prvidx][ptidx] # each trajectory in path
                p.path.append(tr) # store full path
                #print('adding traj: ', tr)
                #
                # sanity check after start
                if iterct > 0:
                    pcurr = p.path[-1].p1
                    pprev = p.path[-2].p2
                    if pprev != pcurr: # check for error
                        print(' prev.p2, curr.p1: ', pprev, pcurr)
                        error(f'path trajectories dont connect! node #: {i}')
            print('testing at iter:',itrct)
            print('p.idxpath: ', p.idxpath)
            #print('p.path: ', p.path)

            c = p.cost()

            print('cost: ', c)
            #x=input('  <cr>')

            if n%10 == 0:
                    print('path ',n, ' cost: ', c) # I'm alive!
            #print('checking path: ', p.idxpath)
            n += 1
            if STOREDATA:
                dfrow = p.idxpath # list of int index pts
                dfrow.append(c)
                dfbf.write(dfrow)
            if c > cmax:
                cmax = c
                pmax.path = p.path
                pmax.idxpath = p.idxpath
                pmax.Tcost = c
                nmax = n
            if c < cmin:
                cmin = c
                pmin.path = p.path
                pmin.idxpath = p.idxpath
                pmin.Tcost = c
                nmin = n
            if n%2000 == 0:
                    print('path ',n, ' cost: ', c) # I'm alive!

            if not LOWMEM:
                path_costs.append(c)
        if not LOWMEM:
            print('quantiling the costs')
            q = [0.25,0.5,0.75,1.0]
            qs = np.quantile(path_costs, q )
            print('Quartile Report:')
            print('  Min:        {:4.1f}'.format( cmin))
            print('    Q1: {:4.2f}  {:4.1f}'.format( q[0], qs[0]))
            print('    Q2: {:4.2f}  {:4.1f}'.format( q[1], qs[1]))
            print('    Q3: {:4.2f}  {:4.1f}'.format( q[2], qs[2]))
            print('    Q4: {:4.2f}  {:4.1f}'.format( q[3], qs[3]))
            print('  Max:        {:4.1f}'.format( cmax ))
        else: # not enough memory for quartiles(!)
            print('Lowest cost path: ', pmin)
            print('path cost: ', cmin)
            print('Highest cost path: ', pmax)
            print('path cost: ', cmax)

        if STOREDATA:
            dfbf.metadata.d['Min Cost']=cmin
            dfbf.metadata.d['Max Cost']=cmax
            if not LOWMEM:
                dfbf.metadata.d['Quartiles']=list(qs)
            dfbf.close()
        pmin.datafile = self.datafile
        #return path object, float
        print('cost set after BF search: ',pmin.costset)
        x = input('... cont')
        return pmin, pmin.Tcost

    def multiHSearch(self,dfile,nsearch):
        ts1 = datetime.datetime.now()
        self.datafile = dfile #keep track of this for adding metadata
        df = dfile
        print('Saving permutations (paths) to: ',df.name)
        itype = str(type(5))
        ftype = str(type(3.1415))
        tps = [itype]*(N*N)      # path point-index sequence
        tps.append(ftype) # the path cost's type
        names = []
        for i in range(N*N):
            names.append('p{:}'.format(i))
        names.append('Cost')
        df.metadata.d['Types'] = tps
        df.metadata.d['Names'] = names
        df.metadata.d['Ncols'] = len(names)
        df.metadata.d['CostType'] = costtype # 'energy' or 'time'
        df.metadata.d['SearchType'] = self.searchtype # 'exhaustive', 'heuristic' etc.
        df.metadata.d['#samples'] = nsearch

        #
        df.open()  # let's open the file (default is for writing)

        # perform and save nsearch heuristic searches
        # from different starting points
        pmin = None
        pmax = None
        cmin = 99999999999
        cmax = 0
        maxTies = 0
        if nsearch > 2*N*N:  # for big enough searches, allocate same # to all start points
            nperstart = nsearch//(N*N)
            MULTI_SEARCH_PER_PT = True
            print(f'searching the {N*N} start points {nperstart} times each.')
        else:
            nperstart = nsearch   # if less, just pick random start points
            MULTI_SEARCH_PER_PT = False
            print(f'searching {nperstart} random start pts, once each.')

        x = input('.....  db pause...')
        for i in range(N*N-1): # go through the start pts
            if MULTI_SEARCH_PER_PT:
                if i%2000==0:
                    print('MH Search iteration: ',i)  #I'm alive
            else:
                pass
            for m in range(nperstart): # do each start pt this many times
                # reset search info
                if not MULTI_SEARCH_PER_PT:
                    # a random start point
                    startPtIdx = random.randint(0,N*N-1)
                else:
                    startPtIdx = i
                print(f'TF: {MULTI_SEARCH_PER_PT} startpt: {startPtIdx} iteration: {m}')
                self.mark[startPtIdx] = False # mark our starting point
                self.Tcost = 0.0
                count = 0

                # do the search
                pself,c = self.heuristicSearch(startPtIdx) #including random tie breakers

                if pself.maxTiesHSearch > maxTies:
                    maxTies = pself.maxTiesHSearch
                datarow = pself.idxpath
                #print('my path: ',datarow)
                datarow.append(c) # last col is cost
                df.write(datarow)
                if c < cmin: # find lowest cost of the runs
                    cmin=c
                    pmin = path(self.grid,self.Cm)
                    pmin.path = pself.path
                    pmin.Tcost = c
                if c > cmax: # find highest cost
                    cmax = c
                    pmax = path(self.grid,self.Cm)
                    pmax.path = self.path
                    pmax.Tcost = c
            if not MULTI_SEARCH_PER_PT:
                break
        df.metadata.d['Min Cost']=cmin
        df.metadata.d['Max Cost']=cmax
        df.metadata.d['Max Ties']=maxTies
        print('min/max path cost: ', cmin,cmax)
        #print('Lowest cost path: ', pmin)
        #print('Highest cost path: ', pmax)
        print('Max # of ties: ',maxTies)

        # save a csv file for tie histogram
        MULTIHISTO = True
        if MULTIHISTO:
            destfolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
            if destfolder[-1] != '/':
                destfolder.append('/')
            tfname = destfolder+'ties_info_'+ df.hashcode+'.csv'
            fp = open (tfname, 'w')
            #print('Distribution of tie choices: (',len(self.path),' points in path)',file=fp)
            #print('\n  tie rank    |  how many times',file=fp)
            #print('----------------------------------------',file=fp)
            #sum = 0
            #medianflag = True
            fmtstring1 = '{:8d} , {:8d} '
            #fmtstring2 = '{:8d}     |     {:8d} << median'
            median=0
            for i,n in enumerate(self.tie_freq):
                median += n
            median /=2
            df.metadata.d['Median Ties']=median
            for i,n in enumerate(self.tie_freq):
                if i>1:
                    if int(n) != 0:
                        logval = np.log10(float(n))
                    else:
                        logval = -1
                    #print(fmtstring1.format(i,int(n),logval),file=fp)
                    print(f"{i:8}   |    {int(n):8}    |    {logval:8.2f}",file=fp)
            fp.close()


        df.close()
        # return path object, float
        return pmin,pmin.Tcost

    def heuristicSearch(self,startptidx):  # path class
        # add to self.path[] one traj at a time greedily
        sr, sc = idx2ij(startptidx)
        crow = sr*N+sc  # starting point in cost matrix
        firstrow = crow # just the row #
        maxTies = 0 # keep track of highest # of tie costs
        # for mh search we can just re-use idxpath and path until done
        self.idxpath = []  # list if index points
        self.path = [] # list of trajectories
        self.Tcost = 0.0 # total path cost
        self.mark = [True for x in range(N*N)]
        if costtype not in ['time','energy']:
            error('unknown cost type: '+costtype)
        #
        #   build the path by greedy algorithm
        #
        while len(self.path) < N*N: # build path up one pt at a time
            # these store all unvisited,valid  branches out of this point/node
            edge_next_tr = []   # next point by trajectory
            edge_costs = []  # cost of branch/traj to next point
            #
            # look at cost of all unmarked,valid branches out of this node
            #
            for ccol in range(N*N):  # ccol is an index
                if self.mark[ccol]:  # only unvisited
                    # look at all unvisited branches leaving current pt
                    br_traj = self.Cm.m[crow][ccol] #branch trajectory from this node
                    #print('\ncrow,ccol:',crow,ccol)
                    #br_traj.aprtr('pull tr from Cm')
                    if br_traj.valid: # don't do self transitions
                        if costtype == 'energy':
                            br_cost = br_traj.e_cost # now pre-computed
                        else:
                            br_cost = br_traj.t_cost
                        #print(' branch cost:',br_cost)
                        edge_next_tr.append(br_traj) # collect all branches out
                        edge_costs.append(br_cost)   # costs of these branches

            # now we have to choose a random branch having min cost
            if len(edge_next_tr)==0:
                error('somethings wrong: i cant find a next node!')
            minCost = min(edge_costs) # will be either a time or energy cost
            epsilon = 0.02*minCost  # within 2% is a tie
            tiebreakerlist = []
            costTmpList = []

            for i,t in enumerate(edge_next_tr): # go through all traj's leaving this pt
                if abs(edge_costs[i]-minCost)<epsilon: # if it's close to min
                    tiebreakerlist.append(t)
                    costTmpList.append(edge_costs[i]) #costs of the tied traj's

            ## debugging junk
            if False and len(tiebreakerlist)>1:
                depth = len(self.idxpath)
                print(depth,': tie: ',len(tiebreakerlist), tiebreakerlist)
                print('costs: ',costTmpList)
                print(' minCost: ', minCost, 'epsil:', epsilon)
                if False:
                    print('current next traj list: ',len(edge_next_tr),'entries')
                    #print(edge_next_tr)
                    print('current cost list:')
                    print(sorted(edge_costs))
                    print('minCost:', minCost)
                    x = input('?...')
                    print('\n\n')
            #####
            #  Collect stats on the tiebreakerlist
            #####
            L = len(tiebreakerlist)
            #error condition checks
            if L < 1:
                error('search.select_next: no next trajs identified yet.')
            if L>maxTies:
                maxTies = L
            if L > MAXTIEHISTO-1:  # pile all bigger ties into last bin
                L = MAXTIEHISTO-1
            self.tie_freq[L] += 1  # count how many ties with each multiplicity L

            # pick a random entry from the ties
            nexttraj = random.choice(tiebreakerlist)
            nextidx = pt2idx(nexttraj.p2) #index of next point (p1 is current pt)

            # Some error checks here
            if not self.mark[nextidx]:
                error('new path-point is marked already')
            if not nexttraj.valid:
                print('crow', crow, 'ccol: ', ccol)
                error('path: invalid new trajectory')
            if crow != firstrow:  # no trajectory has t.p2=startPoint
                pprev = self.path[-1].p2
                pcurr = nexttraj.p1
                if pprev != pcurr: # check for error
                    print('crow/firstrow: ', crow, firstrow)
                    print('adding traj: ', t)
                    print(' prev.p2, curr.p1: ', pprev, pcurr)
                    error('path trajectories dont connect! '+str(crow))

            # OK error checks passed
            self.mark[nextidx] = False  # do not visit this point again
            self.path.append(nexttraj)
            self.idxpath.append(nextidx)
            #print('    adding traj to path: ', self.path[-1])

            # nextidx now become the current point
            crow = nextidx
            self.maxTiesHSearch = maxTies # save this (multi-hsearch will max(max(ties)))
        # don't forget the last point in the idxpath
        self.idxpath.append(pt2idx(self.path[-1].p2))
        print('2D heuristic path search completed!')
        #print('{:} Total path cost ({:}) = {:8.2f}: '.format(self.searchtype,costtype,self.Tcost))
        #print('idxpath: ', self.idxpath)
        #return path object, float
        self.T_cost = self.cost()  # compute and store traj cost
        time.sleep(0.5)
        return self, self.T_cost

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
        if self.datafile is not None:
            hashcode = self.datafile.hashcode
        else:
            hashcode = ''
        fig = plt.figure()
        plt.title('{:}'.format(note))
        plt.xlabel('X     ('+hashcode+')')
        plt.ylabel('\dot{X}')
        plt.grid(True)
        return fig

    def plotOnePath(self,fig):
        # get endpoints for arrows
        x_values = [traj.p1.x for traj in self.path]  # starting values
        y_values = [traj.p1.v for traj in self.path]
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

        # get curves for each arc
        cx, cy = self.compute_curves(-1) #compute trajectory path
        ax.plot(cx,cy,color='blue')
        axlim = 2
        ax.set_xlim([-axlim,axlim])
        ax.set_ylim([-axlim,axlim])

    def plotDone(self):
        plt.show()
        df = self.datafile
        template = f'______________{df.hashcode}.png'
        print('filename template: ',template)
        nroot = input('enter name root: (<enter> to not save) ')
        if len(nroot)>0:
            nroot += '_' # separate the hash
            imgdir = df.folder+'writing/'
            imgname = nroot + df.hashcode
            my_dpi = 200
            plotSave(plt.gcf(), my_dpi, imgdir, imgname)

        else:
            print('plot image NOT saved')




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

