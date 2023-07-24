#!/usr/bin/python3

import socket
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import brl_data.brl_data as bd
import datetime
import sys
import os
import random

from pympler import asizeof

def sizer(str,arg):
    print(' ... click ...')
    #print('Sizeof: ',str, asizeof.asizeof(arg))



def error(msg):
    print('Error: ')
    print(msg)
    quit()

PCNAME = str(socket.gethostname())
N = 4
Npts = N**6
M = Npts

AMAX = 2  #  normalize for now.
DT_TEST = 2.0
DT_START = 1.0  # this needs to be 'smaller' so that Amax can be searched by
                # lengthening dt, but if too small constrain_A can be too slow.

NPC = 20  # number plot points per curve

costtype = 'time'
costtype = 'energy'

# for 1D optimization case
startrow = 3
startcol = 1

def configure(fp=None):
    #print('Pyton path: ', sys.path)
    global costtype, AMAX, DT_TEST, N, M, NPC, startrow, startcol
    if not fp:
        f = open('ctoConfig.txt','r')
    else:
        f = fp

    for line in f:
        parname, *rest = line.split(' ')
        parname = parname.strip()
        v = ' '.join(rest).strip() # parameter value as string

        if parname.startswith('#'):
            next

        if parname =='costtype':
            if v not in ['time','energy']:
                error('path.configure: unknown cost type: '+v)
            costtype = v
        if parname == 'amax':
            AMAX = float(v)
        if parname == 'dt_start': # starting val for constrain_A()
            DT_START = float(v)
        if parname == 'dt_test':
            DT_TEST = dt_test  # fixed dt value used for testing
        if parname == 'N':
            N = int(v)
            M = N**6  # number of grid points
        if parname == 'NPC':  #n time pts in trajectory
            NPC = int(v)
        if parname == 'startrow':
            startrow = int(v)
        if parname == 'startcol':
            startcol = int(v)


#def ij2idx(i,j):
    #error('cant use ij2idx anymore with 6D search')
    #return i*N+j

#def idx2ij(idx):
    #error('cant use idx2ij anymore with 6D search')
    #i = idx//N
    #j = idx-N*i
    #return i,j

#  convert between 6D int coordiates and point index
def getidx(v):
    return v[0]*N**5+v[1]*N**4+v[2]*N**3+v[3]*N**2+v[4]*N+v[5]

def getcoord(idx):
    v = [0,0,0,0,0,0]
    r = idx
    for i in range(6):
        v[5-i] = r%N
        r = (r-v[5-i])//N
    return v

def predict_timing(df, searchtype, systemName, nsamp):
    #
    # predict the search timing
    #
    savedir = ''
    if df is not None:
        savedir = df.folder
        print('set timing data folder to ',df.folder)
    filename = savedir + 'searchTiming.json'
    configKey = str(M) + '-' + searchtype + '-' + systemName
    OK = True
    if os.path.isfile(filename):
        fd = open(filename,'r')
        d = json.load(fd)
        print('n samples:',nsamp)
        try:
            t = d[configKey]
        except:
            OK=False
        if OK:
            print('your predicted search time is:')
            print('search type/PC:',configKey)
            print('rate:          ',d[configKey],'/sec')
            sec = float(d[configKey])*nsamp
            mins = sec/60
            hrs = mins/60
            days = hrs/24
            years = days/365
            print('predicted time: ')
            fmth='{:30} {:>12} {:>12} {:>12} {:>12}'
            print(fmth.format('PC type','mins','hrs','days','years'))
            fmts='{:30} {:12.2f} {:12.2f} {:12.2f} {:12.2f}'
            print(fmts.format(configKey, mins,hrs,days,years))
            if mins>2:
                x=input('OK to continue? ...')
    else:
        OK=False #path is not a file
    if not OK:
        print('no search speed info available for your configuration: ',configKey)
        x=input('OK to continue? ...')

def save_timing(df, searchtype, systemName,rate):
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
    configKey = str(M) + '-' + searchtype + '-' + systemName
    d[configKey] = rate
    fd = open(filename,'w')
    json.dump(d,fd,indent=4)
    return

def plotSave(fig, dpi, imagedir, imagename):
    #fig = plt.gcf()
    #datad = string path to data directory
    if imagedir[-1] != '/':
        imagedir += '/'
    imagepath = imagedir+imagename+'.png'
    fig.savefig(imagepath,dpi=dpi)
    print('your plot is saved to: ',imagepath)

    ####  keep a "log book" if saved
    dim = '6D'
    notes = '{:}, {:}'.format(dim,imagepath)
    now = datetime.datetime.now()
    dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
    logentry = '{:}, plot saved: {:}'.format(dtstring,notes)
    f =open('search_logbook.txt','a')
    print(logentry,file=f)
    f.close()
    ####

class point3D:

    def __init__(self,iv):
        if len(iv) != 6:
            error('point3D: invalid index vector')
        if (N<3):
            error('point3D class: N must be greater than 2')
        self.ix = iv[0]
        self.iy = iv[1]
        self.iz = iv[2]
        self.ixd = iv[3]
        self.iyd = iv[4]
        self.izd = iv[5]
        self.x =     2*iv[0]/(N-1) - 1
        self.y =     2*iv[1]/(N-1) - 1
        self.z =     2*iv[2]/(N-1) - 1
        self.xd = -1*(2*iv[3]/(N-1) - 1)
        self.yd = -1*(2*iv[4]/(N-1) - 1)
        self.zd = -1*(2*iv[5]/(N-1) - 1)
        #self.ivect = [self.ix, self.iy, self.iz, self.ixd, self.iyd, self.izd]
        self.ivect = iv
        if max(self.ivect) > N-1:
            error('point3D: an index is too big for N!')
        self.xvect = [self.x,self.y,self.z,self.xd,self.yd,self.zd]
        self.valid = True

    def randomize(self):
        #
        #  random point locations instead of grid
        #
        self.x =  np.random.uniform(-1,1)
        self.y =  np.random.uniform(-1,1)
        self.z =  np.random.uniform(-1,1)
        self.xd =  np.random.uniform(-1,1)
        self.yd =  np.random.uniform(-1,1)
        self.zd =  np.random.uniform(-1,1)
        self.xvect = [self.x,self.y,self.z,self.xd,self.yd,self.zd]

    def __eq__(self,x):
        for i,s in enumerate(self.ivect):
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
        #print('trajectory3D: p1,p2: ',p1,p2)
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
        dt = DT_START #start with 'fast' dt
        ni = 0
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.05  # if amax too big, slow down
            else:
                break
        dt *= 0.95  # back off prev opt and finetune
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.01
            else:
                break
        am = self.get_Amax(dt)
        #print('Constrain_A: iterations: ', ni, ' dt = {:6.3f} Amax: {:4.2f}'.format(dt,am), ' AMAX:',AMAX)
        #print('             Amax accuracy: {:6.3f}'.format(abs(AMAX - am)))
        self.dt = dt
        self.constrained = True
        self.compute(self.dt)
        return dt

    # also essentially unchanged!
    # 3D
    def timeEvolution(self,ACC_ONLY = False):  #3D
        if not self.computed and not self.constrained:
            error('Cant compute 3D timeEvolution until trajectory is computed and constrained')
        Np = 20-1
        x = []
        v = []
        a = []
        t = []
        for t1 in range (Np):
            t2 = self.dt*t1/Np
            t.append(t2)
            if not ACC_ONLY:
                x.append(self.x(t2))
                v.append(self.xd(t2))
            a.append(self.xdd(t2))
        t.append(self.dt)
        if not ACC_ONLY:
            x.append(self.x(self.dt))
            v.append(self.xd(self.dt))
        a.append(self.xdd(self.dt))
        ##print('A shape1: ',a.shape())
        #
        #  make nicer for plotting and output
        t = np.array(t).T
        if not ACC_ONLY:
            x = np.array(x).T
            v = np.array(v).T
        a = np.array(a).T
        #print('X shape: (tev)',x.shape)
        if ACC_ONLY:
            return a
        else:
            return t,x,v,a

    #
    #   Energy cost: sum of squared acceleration
    def cost_e(self, a):  #3D
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is computed and constrained')
        c = 0.0
        n = len(a)
        for i in range(3):
            for a1 in a:
                c += self.dt*a1[i]*a1[i]/n
        self.e_cost = c
        return c

    def cost_t(self, a):  # need a parameter for efficient calling though not used(!)
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        self.t_cost = self.dt
        return self.dt


    def __repr__(self):
        return str(self.p1) + ' ---> ' + str(self.p2)

def cost_idxp(typestr, Cm, idxpath):   # compute total cost from a list of indeces
    L = len(idxpath)-1
    if L != N**6-1:
        error('not a correct length path!: '+str(L))
    Tc = 0.0
    for i in range(L):
        i1 = idxpath[i]
        i2 = idxpath[i+1]
        ct,ce = Cm.m[i1][i2]
        if typestr == 'energy':
            Tc += ce
        elif typestr == 'time':
            Tc += ct
        else:
            error('cost_idxp: unknown cost type: '+typestr)
    return Tc

class Cm:  # save memory, Cm.m only contains cost pair ct,ce
    def __init__(self,df=None):
        if M>5000:
            error('too many transitions for Cm! ',M*M)

        # this matrix will hold a trajectory for each transition from row to col
        self.m = [[ 0 for x in range(M)] for y in range(M)]
        self.randgrid = False
        if df is not None:
            df.metadata.d['Random Grid'] = False


    def set_GridRandomize(self,df=None):
        self.randgrid = True
        if df is not None:
            df.metadata.d['Random Grid'] = True

        return

    def fill(self):
        print('starting fill...')
        nf = 0
        for i1 in range(M):  # go through all grid points
            for j1 in range(M):
                nf+=1
                if nf%20000==0:
                    pct = i1/M
                    print('fill is {:.0%} done'.format(pct))
                p1 = point3D(getcoord(i1))
                p2 = point3D(getcoord(j1))
                if  self.randgrid:  # select uniform random point location
                    p1.randomize()
                    p2.randomize()
                t = trajectory3D(p1,p2)
                if (p1 == p2):
                    t.valid = False  # eliminate self transitions
                else:
                    #t.compute(DT_TEST)
                    t.constrain_A()
                    a = t.timeEvolution(ACC_ONLY=True)
                    ce = t.cost_e(a)
                    ct = t.cost_t(a)
                    self.m[i1][j1] = (ct,ce)  # new - save data size
        print('done with fill...')

    def __repr__(self):
        res = 'Cm cost matrix:\n'
        for i in range(M):
            for j in range(M):
                res += ' {:5.1f}'.format(self.m[i][j].cost_t)
            res += '\n'
        return res

class search_from_curr_pt:
    def __init__(self,Mark,path):
        self.costtype = costtype # a global from config file
        self.pstartIdx = None  # last known trajectory point (or initial point)
        self.cmin   = 99999999999
        self.minidx = 0
        self.minTrs = []
        self.minidxs = []
        self.found = False
        self.ties = None   # number of ties found in this set of branches
        self.mark = Mark  # array to mark already chosen pts (True == still available)
        self.path = path

    def iterate(self,N,function):
        for ix in range(N):
                for iy in range(N):
                    for iz in range(N):
                        for idx in range(N):
                            for idy in range(N):
                                for idz in range(N):
                                    ivect = [ix,iy,iz,idx,idy,idz]
                                    index = getidx(ivect)
                                    function(index,ivect)

    def find_cmin(self,index,ivect):  # should be called by iterate as first step.
        if self.mark[index]: # index = Cm.m column
            if self.pstartIdx is None:
                error('find_cmin: start point unknown')
            p1idx = self.pstartIdx
            p2idx = index
            if p1idx != p2idx:
                tc = self.eval_cost(p1idx,p2idx)
                if tc < self.cmin:
                    self.cmin   = tc
                    self.minidx = index
                    self.found = True

    def eval_cost(self,i1,i2):
        tc,ec = self.path.Cm.m[i1][i2]
        if self.costtype == 'energy':
            tc = ec #precomputed
        elif self.costtype == 'time':
            tc = tc # precomputed
        else:
            error('search:set_costtype: unknown cost type (3D): '+costtype)
        return tc

    def select_next(self):
        L = len(self.minidxs)

        #error condition checks
        if L != len(self.minidxs):
            brl_error('min trajs and min indexs are different lengths')
        if L < 1:
            error('search.select_next: no next trajs identified yet. did you run find_all_cminTrs?')
        if L == 1:
            self.mark[self.minidx] = False
            return self.minidx, self.minTrs[0]
        if L > MAXTIEHISTO-1:  # pile all bigger ties into last bin
            L = MAXTIEHISTO-1
        self.path.tie_freq[L] += 1  # count how many ties with each multiplicity L
        # pick one at random
        #print('select_next: choosing a random min traj! from ',len(self.minidxs))
        ti = random.choice(self.minidxs) # index of next traj
        tr = trajectory3D(point3D(getcoord(self.pstartIdx)),point3D(getcoord(ti)))
        self.mark[ti] = False
        return ti,tr #chosen next traj index, chosen next traj trajectory

    def find_all_cminTrs(self,index,v): # find a list of all next pts for which cost ~= cmin
        if not self.found:
            error('search.find_all_cminTrs: somethings wrong, need to find_cmin before find_all')
        epsilon = self.cmin * 0.02 # define 'close'
        if self.mark[index]:
            tc = self.eval_cost(self.pstartIdx,index)  # get cost for this branch
            tr = trajectory3D(point3D(getcoord(self.pstartIdx)),point3D(getcoord(index)))
            if abs(tc-self.cmin) < epsilon:
                self.minTrs.append(tr)
                self.minidxs.append(index)
        if len(self.minTrs) > self.path.nmin_max:
            self.path.nmin_max = len(self.minTrs)
        self.ties = len(self.minidxs) #how many ties at this node??

MAXTIEHISTO = 100  # how many bins for the tie frequency histogram

class path3D:
    def __init__(self,Cm):
        self.Cm = Cm      # cost matrix (actually trajectories)
        sizer('path3D.Cm: ',self.Cm)
        self.sr = startrow
        self.sc = startcol
        self.mark = [True for x in range(N**6)]  # true if pt is UNvisited
        sizer('path3D.mark: ',self.mark)
        self.mark[self.sr*N+self.sc] = False # mark our starting point (can be overridden)
        self.Tcost = 0.0
        self.path = []  # the path as a list of trajectories
        self.idxpath = [] # the path as a list of indices (0..N**6)
        self.searchtype = 'none yet'
        self.datafile = None
        self.maxTiesHSearch = -99999999  # most ties when greedy searching
        self.tie_freq = np.zeros(MAXTIEHISTO)  # histogram of how many ties of each length
        self.ties = 0

    def search(self,searchtype,dfile=None,nsamples=1000,profiler=None):
        predict_timing(dfile, searchtype,PCNAME, nsamples)
        #
        #  start timer
        ts1 = datetime.datetime.now()
        #
        #  select the type of search to do
        #
        self.searchtype = searchtype
        #
        #  single heuristic search (greedy nearest neighbor)
        #
        if searchtype.startswith('heur'):
            p, cmin = self.heuristicSearch()
        #
        #  a population of heuristic searches
        #
        elif searchtype.startswith('multi'):
            if dfile is None:
                error('path.search: multi heuristic search requires a dfile')
            p, cmin = self.multiHSearch(dfile,nsamples, profiler=profiler)
        #
        #  a full brute force search over all paths
        #
        elif searchtype.startswith('brute'):
            if dfile is None:
                error('path.search: brute force search requires a dfile')
            p, cmin = self.bruteForce(dfile=dfile,profiler=profiler)
        #
        # a random sample of paths are cost-evaluated
        #
        elif searchtype.startswith('sampling'):
            if dfile is None:
                error('path.search: sampling search requires a dfile')
            p, cmin = self.sampleSearch(dfile=dfile,nsamples=nsamples,profiler=profiler)
        else:
            error('path.search: unknown search type: ', searchtype)
        #report timing
        ts2 = datetime.datetime.now()
        dt = (ts2-ts1).total_seconds()
        print('seconds per {:} paths: {:}'.format(nsamples, float(dt)))
        print('seconds per path: {:}'.format(float(dt)/nsamples))
        save_timing(dfile,searchtype,PCNAME,float(dt)/nsamples)
        return p,cmin

    def sampleSearch(self,dfile=None,nsamples=977,profiler=None):
        print('Testing: sampling search, nsam:',nsamples)
        return self.bruteForce(dfile=dfile,sampling=True,nsamples=nsamples,profiler=profiler)


    def bruteForce(self,dfile=None,sampling=False,nsamples=0,profiler=None): # path class
        if self.searchtype.startswith('none'):
            self.searchtype = 'brute force'
            sampling = False

        if not sampling:
            error('theres NO way I can do full brute force in 6D!!!')

        #n_all_paths = math.factorial(N**6)
        print('Starting {:} search: N={:}'.format(self.searchtype,N))
        if sampling:
            print('   Sampling {:} paths out of <humongous #>'.format(nsamples))

        LOWMEM = False
        if N>3 and not sampling:
            LOWMEM = True  # we just can't store that many paths

        if dfile is not None:
            self.datafile = dfile
            STOREDATA = True   # write all perms and costs to a data file
        else:
            STOREDATA = False

        if STOREDATA:
            df = dfile
            print('Saving permutations (paths) to: ',df.name)
            itype = str(type(5))
            ftype = str(type(3.1415))
            tps = [itype]*(N**6)      # path point seq
            tps.append(ftype) # the path cost's type
            names = []
            for i in range(N**6):
                names.append('{:}'.format(i))
            names.append('Cost')
            df.metadata.d['Ncols'] = len(names)
            df.metadata.d['Types'] = tps
            df.metadata.d['Names'] = names
            df.metadata.d['CostType'] = costtype
            df.metadata.d['SearchType'] = self.searchtype
            df.metadata.d['#samples'] = nsamples
            #
            df.open()  # let's open the file (default is for writing)

        ##1) list all possible paths

        #
        #  Only sampling - full brute force is deleted
        #
        if not sampling:
            error(' we already checked this !!!')

        print('We are generating {:} random paths through {:} nodes'.format(nsamples,N**6))
        phset = set()
        piter = []
        n = 0
        while n < nsamples: # make sure list has no dupes
            if n%10000 == 0:
                print(n,' paths')
            p = list(range(N**6))
            random.shuffle(p) # generate a path as random list of indices
            ph = ''
            for i in range(50):
                ph += str(p[i])[-1] # last digit of idx
            if ph not in phset: # we've found a new pt
                phset.add(ph)
                piter.append(p)
                n+=1 # count adds (faster than len()??)

        sizer('piter: ',piter)

        print('Path enumeration complete (without duplicates):')

        # keep around for history
        #secPerLoop = 0.0003366 # measured on IntelNUC
        #secPerLoop = 0.0008419 # Dell XPS-13

        #2) evaluate their costs
        path_costs = []
        n = -1
        cmin = 99999999999
        cmax = 0
        pmax = path3D(self.Cm)
        sizer('pmax1: ',pmax)
        pmin = path3D(self.Cm)
        sizer('pmax2: ',pmax)
        sizer('pmin: ',pmin)
        for p in piter:  # piter iterates to a series of lists of point indices
            # p is the current path [idx0,idx1,idx2 ...]
            # now get the cost of p
            idxpath = list(p)
            c = cost_idxp(costtype, self.Cm, idxpath)  #what is cost of this path?
            if n%2000 == 0:
                    print('path ',n) # I'm alive!
            if STOREDATA:
                row = idxpath # list of int index pts
                row.append(c)
                df.write(row)
            if c > cmax:
                cmax = c
                #pmax.path = tmpTrajList # list of trajectories
                pmax.idxpath = idxpath
                pmax.Tcost = c
                nmax = n
            if c < cmin:
                cmin = c
                #pmin.path = tmpTrajList
                pmin.idxpath = idxpath
                pmin.Tcost = c
                nmin = n
            #print(' path cost: {:4.2f}'.format(c))
            if not LOWMEM:
                path_costs.append(c)
            if nsamples < 100:  # too much clutter for big searches
                print('{:} path cost: {:12.2f}'.format(n,c))
                sizer('pmin: ',pmin)
        #
        #  we are done with the path set to be evaluated
        #
        print('{:} paths have been evaluated'.format(n))
        print('Lowest cost path: ', pmin)
        print('path cost: ', cmin)
        print('Highest cost path: ', pmax)
        print('path cost: ', cmax)

        if STOREDATA:
            df.metadata.d['Min Cost']=cmin
            df.metadata.d['Max Cost']=cmax
            df.close()
        pmin.datafile = self.datafile
        #return path object, float
        return pmin, pmin.Tcost     # end of bruteforce


    def multiHSearch(self,dfile,nsearch,profiler=None):
        ts1 = datetime.datetime.now()
        self.datafile = dfile #keep track of this for adding metadata
        df = dfile
        print('Saving permutations (paths) to: ',df.name)
        names = []
        for i in range(N**6):
            names.append('{:}'.format(i))
        names.append('Cost')
        itype = str(type(5))
        ftype = str(type(3.1415))
        tps = [itype]*(N**6)      # path point-index sequence
        tps.append(ftype) # the path cost's type (float)
        df.metadata.d['Types'] = tps
        df.metadata.d['Names'] = names
        df.metadata.d['Ncols'] = len(names)
        df.metadata.d['CostType'] = costtype # 'energy' or 'time'
        df.metadata.d['SearchType'] = self.searchtype # 'brute force', 'heuristic' etc.
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
        if nsearch > N**6:  # for big enough searches, allocate same # to all start points
            nperstart = nsearch//N**6
            USESTPT = True
        else:
            nperstart = nsearch   # if less, just pick random start points
            USESTPT = False
        for i in range(N**6): # go through the start pts
            if USESTPT:
                if i%2000==0:
                    print('multiple heuristic searches: ',i)  #I'm alive
            else:
                print('searching starting point:',i)
            # go through the N^2 start points with equal number at each
            startPtIdx = i
            for m in range(nperstart): # do each start pt this many times
                # reset search info
                print('        iteration  ',m,'/',nperstart,'  for starting point',i)
                self.mark = [True for x in range(N**6)]
                count = 0
                if not USESTPT:
                    # a random start point
                    startPtIdx = random.randint(0,N**6-1)

                # don't think we need these acth
                self.mark[startPtIdx] = False # mark our starting point
                self.Tcost = 0.0

                # do the search
                pself,c = self.heuristicSearch3D(startPtIdx) #including random tie breakers
                if maxTies < pself.ties:
                    maxTies = pself.ties
                #
                datarow = pself.idxpath  # ********
                #
                #print('my path: ',datarow)
                datarow.append(c) # last col is cost
                df.write(datarow)
                if c < cmin: # find lowest cost of the runs
                    cmin=c
                    pmin = path3D(self.Cm)
                    pmin.path = pself.path
                    pmin.Tcost = c
                if c > cmax: # find highest cost
                    cmax = c
                    pmax = path3D(self.Cm)
                    pmax.path = self.path
                    pmax.Tcost = c
            if not USESTPT:
                break
        df.metadata.d['Min Cost']=cmin
        df.metadata.d['Max Cost']=cmax
        df.metadata.d['Max Ties']=self.maxTiesHSearch
        print('Lowest cost path: ', pmin)
        print('path cost: ', cmin)
        print('Highest cost path: ', pmax)
        print('path cost: ', cmax)
        print('Max # of ties: ',self.maxTiesHSearch)
        MULTIHISTO = True
        if MULTIHISTO:
            if df.folder[-1] != '/':
                df.folder.append('/')
            fname = df.folder+'ties_info_'+ df.hashcode+'.txt'
            fp = open (fname, 'w')
            print('Distribution of tie choices: (',len(self.path),' points in path)',file=fp)
            print('\n  tie rank    |  how many times',file=fp)
            print('----------------------------------------',file=fp)
            sum = 0
            medianflag = True
            fmtstring1 = '{:8d}     |     {:8d} '
            fmtstring2 = '{:8d}     |     {:8d} << median'
            median=0
            for i,n in enumerate(self.tie_freq):
                median += n
            median /=2
            df.metadata.d['Median Ties']=median
            for i,n in enumerate(self.tie_freq):
                sum += int(n)
                if sum < median and medianflag :
                    fmt = fmtstring1
                elif medianflag:
                    fmt = fmtstring2
                    medianflag = False
                else:
                    fmt = fmtstring1
                print(fmt.format(i,int(n)),file=fp)
            fp.close()

        df.close()
        # return path object, float
        return pmin,pmin.Tcost

    def heuristicSearch3D(self, idx1,profiler=None):
        ADVANCED = True
        BASIC = not ADVANCED  # where' just not going to do BASIC anymore

        # sanity check!!
        if N**6 > 1.0E4:
            error('too big a search!!: '+float(N**6))

        pstartIdx = idx1  # starting point (<N!)
        self.mark = [True for x in range(N**6)]
        count = 0
        self.mark[idx1] = False # mark our starting point


        if ADVANCED:
            self.Tcost = 0.0
            self.path = []
            self.idxpath = [pstartIdx] # path has a start point
            self.nmin_max = 0
            while len(self.path) < N**6 - 1:
                search = search_from_curr_pt(self.mark,self)
                search.pstartIdx = pstartIdx
                search.minTrs=[] #these will get all branches matching cmin cost.
                search.minidxs=[]
                search.iterate(N,search.find_cmin)
                search.iterate(N,search.find_all_cminTrs)
                if search.ties > self.maxTiesHSearch:  # keep track of greatest number of ties along this path
                    self.maxTiesHSearch = search.ties
                nxtidx,nxtTr = search.select_next()   # break a possible tie btwn branches leaving this pt.
                self.path.append(nxtTr)
                self.idxpath.append(nxtidx)
                pstartIdx = nxtidx
                #self.Tcost += search.cmin
            self.Tcost = cost_idxp(costtype,self.Cm, self.idxpath)  # compute total cost of the path

            REPORTHISTO = False
            #if REPORTHISTO:
                #print('Distribution of tie choices: (',len(self.path),' points in path)')
                #print('\n  tie rank    |  how many times')
                #print('----------------------------------------')
                #sum = 0
                #medianflag = True
                #fmtstring1 = '{:8d}     |     {:8d} '
                #fmtstring2 = '{:8d}     |     {:8d} << median'
                #median=0
                #for i,n in enumerate(self.tie_freq):
                    #median += n
                #median /=2

                #for i,n in enumerate(self.tie_freq):
                    #sum += int(n)
                    #if sum < median and medianflag :
                        #fmt = fmtstring1
                    #elif medianflag:
                        #fmt = fmtstring2
                        #medianflag = False
                    #else:
                        #fmt = fmtstring1
                    #print(fmt.format(i,int(n)))
                #print('')
                #print('\n\n       The longest set of min-cost next points was: {:} points\n\n'.format(self.nmin_max))
            print('heuristic path search completed!')
            print('{:} Total path cost ({:}) = {:8.2f}: '.format(self.searchtype,costtype,self.Tcost))
            # return path object, float
            return self, self.Tcost

        if BASIC:  # we're not doing basic anymore
            error("we're not doing BASIC mode anymore")

    def check(self): # 3D
        if len(self.path) != N**6-1:
            error('wrong path length '+str(len(self.path)))
        i=0
        for t in self.path:
            if not t.valid:
                error('path contains invalid traj: '+str(i))
            if i > 0:
                if self.path[i-1].p2 != t.p1:
                    error('discontinuous path: '+str(i-1)+' '+str(i))
            i+=1


    def compute_curves(self,idx): #3D
        curvepts_x = []
        for i,tr in enumerate(self.path):
            if i == idx or idx < 0:
                if i == len(self.path):
                    break
                if not tr.valid:
                    error('compute_curves()3D: I should not have found an invalid trajectory in path: '+str(i))
                if tr is None:
                    error('null traj: '+ str(i) + str(tr))
                if not tr.computed and not tr.constrained:
                    error('Cant plot until trajectory is constrained: '+ str(i))
                dt = tr.dt
                for i in range(NPC):
                    t = dt*i/NPC
                    curvepts_x.append(tr.x(t))
                curvepts_x.append(tr.x(dt))
        x = np.array(curvepts_x).T
        return x


    def save(self,fname): # save 3D trajectory for animation and plotting
        #
        #   this is a new datafile just for visualization
        #
        df = bd.datafile(fname, 'BH', 'simulation')
        df.set_folders('','') # default local foldesr
        df = self.datafile
        x = self.compute_curves(-1) # save all trajectories
        print('path.save: x points:     ',len(self.path))
        print('path.save: x curves dims:', x.shape)
        col_names = ['n','X','Y','Z']
        int_type = str(type(5))    # these have to be strings b/c json can't serialize types(!)
        float_type = str(type(3.14159))
        col_types = [ int_type, float_type, float_type, float_type]
        col_comments = ['','','','']
        df.set_metadata(col_names, col_types, col_comments)
        df.metadata.d['AMAX']=AMAX
        df.metadata.d['N'] = N
        df.metadata.d['costtype'] = costtype
        df.metadata.d['Searchtype'] = self.searchtype
        df.metadata.d['data purpose'] = 'Visualization/Animation'
        if df.validate():
            print('Datafile is properly set up with valid metadata')
        else:
            print('somethings wrong with datafile or metadata setup')

        df.open()  # let's open the file (default is for writing)

        ##  Now lets write out data
        r,c = np.shape(x)
        for i in range(c):
            row = [i,x[0][i],x[1][i],x[2][i]]
            df.write(row)

        df.close()   # all done


def main():
    import pickle
    import os

    epsilon = 0.02

    configure()

    assert N==3

    print('Commencing tests: 6D, N=',N)

    print('index <--> coordinates test')

    v=[0,0,0,0,0]  # 5 random 6D test coordinate vectors

    for i in range(len(v)):
        v[i] = [0]*6
        for j in range(6):
            v[i][j] = random.choice(range(N))
        print('v: [',i,']:',v[i])

    for i in range(3):
            x = getcoord(getidx(v[i]))
            assert x==v[i]

    print('index <--> coordinates test: PASSED')


    print('\n\n   Cm tests: ')

    c1 = Cm()

    SKIPFILL = False  # save time to focus on later tests

    if not SKIPFILL:
        pname = 'c1Costs.pickle'
        if os.path.exists(pname):
            print('loading precomputed cost matrix   ...')
            f = open(pname, 'rb')
            c1 = pickle.load(f)
            print(' loading completed')
        else:
            print('no stored data: computing cost matrix')
            c1.fill()
            f = open(pname,'wb')
            pickle.dump(c1,f)

    r,c = np.shape(c1.m)

    assert r==c
    assert r == N**6

    p1 = point3D(getcoord(1234))
    p2 = point3D(getcoord(3241))
    tr12 = trajectory3D(p1,p2)
    tt3d = type(tr12)

    if not SKIPFILL:
        # [10][10] is just a "random" traj
        assert type(c1.m[10][10]) == tt3d
        assert not c1.m[10][10].valid  # self-self transitions invalid.
        assert c1.m[10][11].valid
    else:
        print('\n                 fill test skipped!!\n')
        c1.m[10][10] = tr12

    print('   Cm tests:  PASSED')

    print('\n\n   Amax tests')
    print('\n\n If this test fails, try a smaller value of DT_START (in ctoConfig.txt)\n')

    print('DT_START: ',DT_START)
    print('AMAX: ',AMAX)
    for i in range(10):
        p1 = point3D(getcoord(random.randint(0,N**6-1)))
        p2 = point3D(getcoord(random.randint(0,N**6-1)))
        tr12.p1 = p1
        tr12.p2 = p2
        tr12.constrain_A()
        print('Amax test: tr12:',tr12)
        a = tr12.timeEvolution(ACC_ONLY=True)
        amax_min = 9999999999
        for i in range(3):
            #print('amax[i]: ',max(a[i]), min(a[i]))
            # is acceleration properly constrained?
            mx = max(abs(max(a[i])), abs(min(a[i])))
            amax_diff = abs(mx-AMAX)/AMAX
            if amax_diff < amax_min:
                amax_min = amax_diff
        print('             amax_min: ',amax_min) # how close peak acc is to AMAX
        # test that amax is close to AMAX
        assert amax_min < epsilon

    print('\n\n   Amax tests:   PASSED')


    print('\n\n    Cost tests')
    a = tr12.timeEvolution(ACC_ONLY=True)
    tr12.cost_e(a)
    tr12.cost_t(a)
    assert tr12.e_cost > 0
    assert tr12.t_cost > 0


    print('cost tests:     PASSED')

    #
    #  set up a df for search tests
    #
    df = bd.datafile('tests','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')
    df.metadata.d['ResearchQuestion'] = 'debug and test'


    print('\n\nsearch test 1:   sampling random paths')
    searchtype = 'sampling search'
    bpath = path3D(c1)
    bpath.search(searchtype,dfile=df,nsamples=10000)


    print('\n\n            ALL tests:     PASSED')

if __name__ ==  '__main__':
    print('main starting:')
    main()

