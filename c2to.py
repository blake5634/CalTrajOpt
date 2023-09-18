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

# all the following will be overridden by searchOpt.py

N = 4
Npts = N**6  # will be overridden by searchOpt SPACE flag
M = Npts
AMAX = 2  #  normalize for now.
DT_TEST = 2.0
DT_START = 1.5  # this needs to be 'smaller' so that Amax can be searched by
                # lengthening dt, but if too small constrain_A can be too slow.
NPC = 20  # number plot points per curve
costtype = 'time'
costtype = 'energy'
# for 1D optimization case
startrow = 3
startcol = 1


#
#  Scale dims for each axis to match a real robot
#  examples:
#  [-1,1] means keep axis at the default normalized values
#  [-30,25] means axis range should be -30-->25
#
# defaults: (no scaling)
Scale2D = [[-1,1], [-1,1]]
Scale6D = [[-1,1], [-1,1], [-1,1], [-1,1], [-1,1], [-1,1]]

MAXTIEHISTO = 80  # how many bins for the tie frequency histogram

def configure(fp=None):
    #print('Pyton path: ', sys.path)
    global costtype, AMAX, DT_TEST, N, M, NPC, startrow, startcol, gridType
    gridType = 'rectangular'
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
                error('configure(): unknown cost type: '+v)
            costtype = v
        if parname =='gridType':
            if v not in ['random','rectangular']:
                error('configure(): unknown gridType: '+v)
            gridType = v
        if parname == 'SPACE':
            SPACE = v
            if v not in(['2D','6D']):
                error('configure(): SPACE must be "2D" or "6D" not ',v)
        if parname == 'amax':
            AMAX = float(v)
        if parname == 'dt_start': # starting val for constrain_A()
            DT_START = float(v)
        if parname == 'dt_test':
            DT_TEST = dt_test  # fixed dt value used for testing
        if parname == 'N':
            N = int(v)
        if parname == 'NPC':  #n time pts in trajectory
            NPC = int(v)
        if parname == 'startrow':
            startrow = int(v)
        if parname == 'startcol':
            startcol = int(v)

    if SPACE == '2D':  # do this after all options read
        Npts = N*N
    else:
        Npts = N**6

#  'i,j' 'r,c', 'row,col' are used kind of interchangeably
#   to indicate point in rectangular grid ... sorry!
def ij2idx(i,j):   #  0-Npts-1
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
    dim = '2D'  # change if on multiOpt branch
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
        self.randgrid = False
        self.fromFile = False
        self.gr = []
        for i in range(self.N):
            row = []
            for j in range(self.N):
                row.append(point2D(i,j))
            self.gr.append(row)

    def scale2D(self):
        for i in range(self.N):
            for j in range(self.N):
                self.gr[i][j].scale2D()

    #2D class grid2D
    def savePoints2D(self,df):   # save points generated into a file
        if not self.randgrid:
            error('savePoints2D: Do not save points if they are not random')
        # each entry in Cm is a trajectory2D  we need only store p1 and p2
        it =str(type(5))
        ft =str(type(3.14159))
        df.metadata.d['Types'] = [it,it,ft,ft]
        df.metadata.d['Names'] = ['i1','j1','x1','v1']
        df.metadata.d['Ncols'] = len(df.metadata.d['Names'])
        df.metadata.d['Space'] = '2D'
        df.metadata.d['Grid dim'] = N
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

    #class grid2D
    def readPoints2D(self,df):  # read randomized points from a file
        if not self.randgrid:
            error('readPoints2D: Do not read in points if they are not random')
        df.open('r')
        # read in pair of points for each Cm.m entry and store trajectory.
        #  random points are U(-1,1).
        ptindex = 0
        for row in df.reader:

            i = int(row[0])  # 0 -- Npts-1
            j = int(row[1])

            print(f"I'm creating random point at {i}, {j}")
            p1 = point2D(i,j)
            self.gr[i][j]=p1
            ptindex += 1
        df.close()
        self.fromFile = True
        # return the file name from which the pts were read for the record
        return df.hashcode

#  convert between 6D int coordinates and point index
def getidx6D(v):
    return v[0]*N**5+v[1]*N**4+v[2]*N**3+v[3]*N**2+v[4]*N+v[5]

def getcoord6D(idx):
    v = [0,0,0,0,0,0]
    r = idx
    for i in range(6):
        v[5-i] = r%N
        r = (r-v[5-i])//N
    return v

def setupPoints6D():
    global pts, Cm
    Cm = Cm6D() # initially all zeros
    pts = []
    if gridType == 'random':
        print(f'TEST:  generating random grid: {Npts} rows')
        for i in range(Npts):
            newpt = point6D(getcoord6D(i))
            newpt.randomize()
            print('r',end='')
            pts.append(newpt)
    else:
        for i in range(Npts):
            newpt = point6D(getcoord6D(i))
            print('.',end='')
            pts.append(newpt)
    print('')
    print(f'setupPoints6D generated {len(pts)} points')
    return pts

def savePoints6D(df):   # save points generated (esp random) into a file
    # each entry in Cm is a trajectory2D  we need only store p1 and p2
    it =str(type(5))
    ft =str(type(3.14159))
    df.metadata.d['Types'] = [it]*6 + [ft]*6
    # Names =    ['i1','j1','x1','v1']
    df.metadata.d['Names'] = [f'c{i}' for i in range(6)]+[f'x{i}' for i in range(6)]
    df.metadata.d['Ncols'] = len(df.metadata.d['Names'])
    if df.metadata.d['Ncols'] != 12:
        error('somethings wrong in savePoints6D')
    df.metadata.d['Space'] = '6D'
    df.descript_str = 'randomGridPointSet' # standardize for point storage filename fields
    df.metadata.d['Research Question'] = 'RandomGridPointSet' # standardize for RQ
    df.open('w')
    # write out key info for each point.
    for idx in range(Npts):
        p1 = pts[idx]  # precomputed points
        #row = [i,j, p1.x, p1.v ]
        row = []
        for i2 in getcoord6D(idx):
            row.append(i2)
        row = p1.ivect
        row += p1.xvect
        df.write(row)
    df.close()

def readPoints6D(df):  # read randomized points from a file
    global pts
    df.open('r')
    pts = []  # place for the data
    # read in a point for each entry store in pts global
    ptindex = 0
    for row in df.reader:
        siv = row[0:6]
        sxv = row[6:]
        if len(siv) != 6 or len(sxv) != 6:
            error('something wrong line 226')
        iv = []
        xv = []
        for s in siv:   # convert from strings returned by reader
            iv.append(int(s))
        if ptindex != getidx6D(iv):
            print(' ... somethings wrong line 297')
        p1 = point6D(iv)
        print(f"I'm creating 6D point at {iv}")
        p1.x = xv[0]  # transfer the random data to this pt.
        p1.y = xv[1]
        p1.z = xv[2]
        p1.xd = xv[3]
        p1.yd = xv[4]
        p1.zd = xv[5]
        p1.xvect = xv
        pts.append(p1)
        ptindex += 1
    df.close()
    # return the file name from which the pts were read for the record
    return df.hashcode

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



class point6D:
    def __init__(self,iv):
        if len(iv) != 6:
            error('point6D: invalid index vector')
        if (N<3):
            error('point6D class: N must be greater than 2')
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
            error('point6D: an index is too big for N!')
        self.xvect = [self.x,self.y,self.z,self.xd,self.yd,self.zd]
        self.valid = True

    def randomize(self):
        #
        #  random point locations instead of grid
        #
        self.x  =  np.random.uniform(-1,1)
        self.y  =  np.random.uniform(-1,1)
        self.z  =  np.random.uniform(-1,1)
        self.xd =  np.random.uniform(-1,1)
        self.yd =  np.random.uniform(-1,1)
        self.zd =  np.random.uniform(-1,1)
        self.xvect = [self.x,self.y,self.z,self.xd,self.yd,self.zd]


    def scale6D(self):
        vr = []
        for i in range(6):
            #generate Xform to scale [-1,1] to desired range
            # y = a*(x+x0) + b
            srange = Scale6D[i][1] - Scale6D[i][0]  # axis 1
            a = srange/2
            b = Scale6D[i][0]
            vr.append(a*float(self.xvect[i]+1.0)+b)
        self.xvect = vr

        self.x  = vr[0]
        self.y  = vr[1]
        self.z  = vr[2]
        self.xd = vr[3]
        self.yd = vr[4]
        self.zd = vr[5]

    def __eq__(self, x):
        for i,sxv in enumerate(self.xvect):
            xxv = x.xvect[i]
            if abs(sxv-xxv) > 0.01:
                return False
        return True

    #6D
    def __repr__(self):
        txt = '['
        for x in self.xvect:
            txt += f' {x:5.1f}'
        txt += ']'
        return txt

class point2D:
    def __init__(self,i,j):  # index the spatial grid 0--N-1, 0--N-1
        if i > N or j>N or i<=-1 or j<=-1:
            error('point2D i,j is too big: '+str(i)+' '+str(j)+ ' /'+str(N))
        self.row = i
        self.col = j
        self.x =     2*j/(N-1) - 1   # these subject to change with .scale2D()
        self.v = -1*(2*i/(N-1) - 1)
        if self.x < -1.05 or self.x > 1.05:
            msg = f"point2D: I'm creating bogus point coordinates: {i} {j} {self.x} {self.v}"
            error(msg)
        self.tr = None

    def randomize(self):
        self.x = np.random.uniform(-1,1)
        self.v = np.random.uniform(-1,1)
        return

    def scale2D(self):
        vr = [self.x, self.v]
        v = []
        for i,val in enumerate(vr):
            #generate Xform to scale [-1,1] to desired range
            # y = a*(x+x0) + b
            a = (Scale2D[i][1] - Scale2D[i][0])/2.0
            b = Scale2D[i][0]
            v.append(a*float(val+1.0)+b)
        self.x  = v[0]
        self.v  = v[1]
        return

    def __eq__(self,x):
        #return self.row == x.row and self.col == x.col
        epsilon = 0.02
        ex = abs(self.x-x.x) < epsilon # range normalized to -1,1
        ev = abs(self.v-x.v) < epsilon
        return ex and ev

    def __repr__(self):
        return ' ({:4.2f}, {:4.2f})'.format(self.x, self.v)
        #return ' ({:}, {:})'.format(self.row, self.col)
        st = ' ('
        for s in self.ivect:
            st += '{:3d}, '.format(s)
        st += ')'
        return st

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
        for a1 in a:
            c += self.dt*a1*a1
        self.e_cost = c/len(a)
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is computed and constrained')
        self.t_cost = self.dt
        return self.dt

    def __repr__(self):
        return str(self.p1) + ' ---> ' + str(self.p2)


class trajectory6D:
    global Cm

    def __init__(self,p1,p2):
        #print('trajectory6D: p1,p2: ',p1,p2)
        self.p1 = p1
        self.p2 = p2
        self.a0 = [0,0,0]
        self.a1 = [0,0,0]
        self.a2 = [0,0,0]
        self.a3 = [0,0,0]
        self.computed = False
        self.constrained = False
        self.valid = True
        self.e_cost = None
        self.t_cost = None

    #6D
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

    #6D
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

    #6D
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
    # it seems constrain_A can be identical for 6D and 1D(!)
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
                dt *= 1.10  # if amax too big, slow down
            else:
                break
        dt *= 0.9  # back off prev opt and finetune
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.04
            else:
                break
        am = self.get_Amax(dt)
        #print('Constrain_A: iterations: ', ni, ' dt = {:6.3f} Amax: {:4.2f}'.format(dt,am), ' AMAX:',AMAX)
        #print('             Amax accuracy: {:6.3f}'.format(abs(AMAX - am)))
        self.dt = dt
        self.constrained = True
        self.compute(self.dt) # because we changed dt
        return dt

    #6D
    def timeEvolution(self,ACC_ONLY = False):  #6D & 6D
        if not self.computed and not self.constrained:
            error('Cant compute 6D timeEvolution until trajectory is computed and constrained')
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
    #6D
    def getCosts(self,Cm):
        idx1 = getidx6D(self.p1.ivect)
        idx2 = getidx6D(self.p2.ivect)
        if Cm.m[idx1][idx2] != 0: # if we've already computed
            ct,ce = Cm.m[idx1][idx2]
        else:
            self.constrain_A()
            a = self.timeEvolution(ACC_ONLY=True)
            ce = self.cost_e(a)
            ct = self.cost_t(a)
            Cm.m[idx1][idx2] = (ct,ce) # store for potential re-use
        self.t_cost = ct
        self.e_cost = ce
        return ct,ce
    #
    #   Energy cost: sum of squared acceleration
    def cost_e(self, a):  #6D
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

#6D
def getCosts_ij(idx1,idx2,Cm):
    if Cm.m[idx1][idx2] != 0: # if we've already computed
        ct,ce = Cm.m[idx1][idx2]
    else:
        tr = trajectory6D(pts[idx1],pts[idx2])
        tr.constrain_A()
        a = tr.timeEvolution(ACC_ONLY=True)
        ce = tr.cost_e(a)
        ct = tr.cost_t(a)
        Cm.m[idx1][idx2] = (ct,ce) # store for potential re-use
    return ct,ce

class Cm2D: # matrix full of trajectory2D objs
    def __init__(self,df=None):
        self.m = [[ 0 for x in range(Npts)] for y in range(Npts)]
        # normally the points are in a regular 6 dim grid.  if Randgrid
        # is true, they will be converted into uniform([-1,1)) in all coordinates
        self.randgrid = False
        if df is not None:
            df.metadata.d['Random Grid'] = False
    #2D
    def set_GridRandomize(self,df=None):
        self.randgrid = True
        if df is not None:
            df.metadata.d['Random Grid'] = True
        return
    #2D
    def fill(self, grid):
        nselftr = 0
        print('starting fill...{:}x{:}'.format(N,N))
        grid.randgrid = False # we do the setting of this and the randomization itself HERE
        if grid.fromFile:  # we have read grid already from file and it must be random
            grid.randgrid=True
        elif self.randgrid:  # but if haven't read a file, (gen from scratch)
            grid.randgrid=True
            # go through the grid and randomize their x,v
            for r in range(N):  # for all 2nd points
                for c in range(N):
                    ptmp = grid.gr[r][c]
                    ptmp.randomize()
                    grid.gr[r][c] = ptmp
        else: # no file and not random: non random grid
            pass # grid should already have x,v coordinates
        # now get a tr (and costs) for all possible transitions (row x col)
        for i1 in range(Npts):
            for j1 in range(Npts):
                r1,c1 = idx2ij(i1)
                r2,c2 = idx2ij(j1)
                p1 = grid.gr[r1][c1]
                p2 = grid.gr[r2][c2]
                # we now have p1 and p2
                t = trajectory2D(p1,p2)
                t.valid = True
                if (p1 == p2): # robust to random pts
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

def cost_idxp6D(typestr, idxpath, Cm):   # compute total cost from a list of indices
    L = len(idxpath)-1
    if L != Npts-1:
        error('cost_idxp6D: not a correct length path!: '+str(L))
    Tc = 0.0
    for i in range(L):
        i1 = idxpath[i]
        i2 = idxpath[i+1]
        if i1==i2:
            error('path repeats a node index')
        #ct,ce = Cm.m[i1][i2]
        tr = trajectory6D(pts[i1],pts[i2])
        tr.getCosts(Cm)
        ce = tr.e_cost
        ct = tr.t_cost
        if typestr == 'energy':
            Tc += ce
        elif typestr == 'time':
            Tc += ct
        else:
            error('cost_idxp: unknown cost type: '+typestr)
    return Tc

class Cm6D:  # save memory, Cm.m only contains cost pair ct,ce
    def __init__(self,df=None):
        if M>5000:
            error('too many transitions for Cm! ',M*M)

        # this matrix will hold a trajectory for each transition from row to col
        self.m = [[ 0 for x in range(M)] for y in range(M)]
        # normally the points are in a regular 6 dim grid.  if Randgrid
        # is true, they will be converted into uniform([-1,1)) in all coordinates
        self.randgrid = False
        if df is not None:
            df.metadata.d['Random Grid'] = False

    #6D
    def set_GridRandomize(self,df=None):
        self.randgrid = True
        if df is not None:
            df.metadata.d['Random Grid'] = True
        return

    #6D
    def fill(self):
        print('starting fill...')
        error('not using fill anymore')
        if self.randgrid:  # need to store a point for each col
            #****************
            #*
            #*    Have to rethink this!!!
            #*********************************
            colpoints = [] # store the randomized point for each col (just once)
            for i2 in range(N):  # for all 2nd points
                for j2 in range(N):
                    p2 = grid.gr[i2][j2]  # i2,j2 have same end point p2
                    p2.randomize()   #only randimize ONCE per row
                    colpoints.append(p2)

        nf = 0
        for i1 in range(M):  # go through all grid points
            for j1 in range(M):
                nf+=1
                if nf%2000==0:
                    pct = i1/M
                    print('fill is {:.0%} done'.format(pct))
                p1 = point6D(getcoord6D(i1))
                p2 = point6D(getcoord6D(j1))
                if  self.randgrid:  # select uniform random point location
                    p1.randomize()  # this is broken and we're not using .fill anymore! (8/2)
                    p2.randomize()
                t = trajectory6D(p1,p2)
                if (p1 == p2):
                    t.valid = False  # eliminate self transitions
                else:
                    #
                    #   save the cost in Cm (old: full trajectory in Cm
                    ##       was big mem hog.)
                    t.getCosts(self) # since we are already in class Cm!
                    #t.constrain_A()
                    #a = t.timeEvolution(ACC_ONLY=True)
                    #ce = t.cost_e(a)
                    #ct = t.cost_t(a)
                    if (p1 == p2):
                        t.valid = False  # eliminate self transitions
                        self.m[i1][j1] = (0,0)
                    else:
                        self.m[i1][j1] = (ct,ce)
        print('done with fill...')
        print('# of self-state (invalid) transitions: ',nselftr, ' (expected Npts): ', Npts)

    def __repr__(self):
        res = 'Cm cost matrix:\n'
        if Npts > 100:
            error('Cm is too big to print out!')
        else:
            for i in range(Npts):
                for j in range(Npts):
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

class path2D:
    def __init__(self,grid,Cm):
        self.Cm = Cm      # cost matrix (actually trajectories)
        self.grid = grid
        self.mark = [True for x in range(Npts)]
        #self.mark[self.sr*N+self.sc] = False # mark our starting point
        self.Tcost = 0.0
        self.path = []  # the path as a list of trajectories
        self.idxpath = [] # the path as a list of indeces (0..Npts)
        self.searchtype = 'none yet'
        self.datafile = None
        #collect tie stats on this path
        self.maxTiesHSearch = -99999999  # most ties when greedy searching
        self.tie_freq = np.zeros(MAXTIEHISTO)  # histogram of how many ties of

    #2D
    def search(self,searchtype,dfile=None,nsamples=1000):
        self.searchtype = searchtype
        print(f'testing: {searchtype}, {dfile.name}, {nsamples}')
        predict_timing(dfile, searchtype, PCNAME, nsamples)
        #
        #  start timer
        ts1 = datetime.datetime.now()
        #
        #  select the type of search to do
        #
        if searchtype.startswith('heur'):
            startidx = 1
            p, cmin = self.heuristicSearch(startidx)
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
        perfdt = (ts2-ts1).total_seconds()
        print('seconds per {:} paths: {:}'.format(nsamples, float(perfdt)))
        print('seconds per path: {:} (includes plot viewing time!)'.format(float(perfdt)/nsamples))
        save_timing(dfile,searchtype,PCNAME,float(perfdt)/nsamples)
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

    #2D
    def sampleSearch(self,dfile=None,nsamples=977):
        return self.bruteForce(dfile=dfile,sampling=True,nsamples=nsamples)

    #2D
    def bruteForce(self,dfile=None,sampling=False,nsamples=0): # path class
        if self.searchtype.startswith('none'):
            self.searchtype = 'exhaustive'
        n_all_paths = math.factorial(Npts)
        print('Starting {:} search: N={:}'.format(self.searchtype,N))
        if sampling:
            print(f'   Sampling {nsamples} paths out of {float(n_all_paths):12.5e}')
        else:
            print(f'   Exhaustive search of {float(n_all_paths):12.5e} paths.')

        x = input(' ... ready?  ...')

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
            tps = [itype]*(Npts)      # path point seq
            tps.append(ftype) # the path cost's type
            names = []
            for i in range(Npts):
                names.append(f'p{i}')
            names.append('Cost')
            dfbf.metadata.d['Ncols'] = len(names)
            dfbf.metadata.d['Space'] = '2D'
            dfbf.metadata.d['Types'] = tps
            dfbf.metadata.d['Names'] = names
            dfbf.metadata.d['CostType'] = costtype
            dfbf.metadata.d['SearchType'] = self.searchtype
            #dfbf.metadata.d['#samples'] = nsamples  SHOULD BE AUTO

            #
            dfbf.open()  # let's open the file (default is for writing)

        ##1) list all possible paths
        if not sampling:
            print('We are about to find all paths through ',Npts,' nodes')
            x=input('ready?..')
            piter = itt.permutations(range(Npts),Npts) # not a list!
        else:
            print('We are generating {:} random paths through {:} nodes'.format(nsamples,Npts))
            piter = []
            phashset = set() # just for detection dupes
            while len(phashset) < nsamples:
                p = list(range(Npts))
                random.shuffle(p) # generate a path as random list of indices
                pthash = ''
                for j in p:
                    pthash +='{:3d}'.format(j) # 'hash' the list
                if pthash not in phashset:
                    piter.append(p) # the actual path
                    phashset.add(pthash) # robust to random duplicates

        print('Path enumeration complete:')

        # keep around for history
        #secPerLoop = 0.0003366 # measured on IntelNUC
        #secPerLoop = 0.0008419 # Dell XPS-13

        #2) evaluate their costs
        path_costs = []
        n = -1
        cmin = 99999999999
        cmax = 0
        pmax = path2D(self.grid,self.Cm)   # storage for int results
        pmin = path2D(self.grid, self.Cm)

        itrct = -1
        for p in piter:  # piter returns a whole "list" of point indices each time
            itrct += 1
            idxpath = list(p)
            p = path2D(self.grid, self.Cm) # temp variable
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
                # valid flag check:
                if prvidx != ptidx and not tr.valid:
                    print('Ive got a wierd .valid flag (should be True)')
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
        return pmin, pmin.Tcost

    #2D
    def multiHSearch(self,dfile,nsearch):
        ts1 = datetime.datetime.now()
        self.datafile = dfile #keep track of this for adding metadata
        df = dfile
        print('Saving permutations (paths) to: ',df.name)

        # set up the output file metadata
        itype = str(type(5))
        ftype = str(type(3.1415))
        tps = [itype]*(Npts)      # path point-index sequence
        tps.append(ftype) # the path cost's type
        names = []
        for i in range(Npts):
            names.append('p{:}'.format(i))
        names.append('Cost')
        df.metadata.d['Types'] = tps
        df.metadata.d['Names'] = names
        df.metadata.d['Ncols'] = len(names)
        df.metadata.d['CostType'] = costtype # 'energy' or 'time'
        df.metadata.d['Space'] = '2D'
        df.metadata.d['Grid dim'] = N
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

        if nsearch > 2*Npts:  # for big enough searches, allocate same # to all start points
            nperstart = nsearch//Npts
            MULTI_SEARCH_PER_PT = True
            loopcount = Npts
        else:
            nperstart = 1   # if less, just pick random start points
            MULTI_SEARCH_PER_PT = False
            loopcount = nsearch
        print(f'starting {nsearch} 2Diterations with loopcount {loopcount} and nerstart {nperstart}')
        x = input (' ... testpause ...')
        for i in range(loopcount): # go through the start pts
            for m in range(nperstart): # do each start pt this many times
                # reset search info
                if not MULTI_SEARCH_PER_PT:
                    # a random start point
                    startPtIdx = random.randint(0,Npts-1) #randint includes the 2nd val!
                else:
                    startPtIdx = i
                self.Tcost = 0.0
                count = 0

                print(f'       2D iteration  {m}/{nperstart}  for starting point {i+1}/{loopcount}')
                # do the search
                pself,c = self.heuristicSearch(startPtIdx) #including random tie breakers

                if pself.maxTiesHSearch > maxTies:
                    maxTies = pself.maxTiesHSearch
                #
                datarow = pself.idxpath
                #
                #print('my path: ',datarow)
                datarow.append(c) # last col is cost
                df.write(datarow)
                if c < cmin: # find lowest cost of the runs
                    cmin=c
                    pmin = path2D(self.grid,self.Cm)
                    pmin.path = pself.path
                    pmin.Tcost = c
                if c > cmax: # find highest cost
                    cmax = c
                    pmax = path2D(self.grid,self.Cm)
                    pmax.path = self.path
                    pmax.Tcost = c
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

    #2D
    def heuristicSearch(self,startptidx):  # path class
        # add to self.path[] one traj at a time greedily

        #crow indicates the point we are incrementally searching from
        #     possible branches from crow are the cols, ccol (if not yet added
        #     and not self transitions)
        crow = startptidx  # first crow is here
        firstrow = crow # just the row number

        maxTies = 0 # keep track of highest # of tie costs
        self.idxpath = [startptidx]  # list if index points
        self.path = [] # list of trajectories
        self.Tcost = 0.0 # total path cost
        self.mark = [True for x in range(Npts)]
        self.mark[startptidx] = False  # our first pt must be marked
        if costtype not in ['time','energy']:
            error('unknown cost type: '+costtype)
        #
        #   build the path by greedy algorithm
        #
        print('   >>> New hsearch')
        #  correct length of path is Npts-1 trajectories
        while len(self.path) < Npts-1: # build path up one Traj at a time
            # these store all unvisited,valid  branches out of this point/node
            edge_next_tr = []   # trajectory to the next pt
            edge_costs = []  # cost of branch/traj to the next point
            #
            # look at cost of all unmarked,valid branches out of this node
            #
            for ccol in range(Npts):  # ccol is an index
                if self.mark[ccol]:  # only unvisited
                    # look at all unvisited branches leaving current pt
                    br_traj = self.Cm.m[crow][ccol] #branch trajectory from this node
                    #print('\ncrow,ccol:',crow,ccol)
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
            nexttraj =  random.choice(tiebreakerlist)
            nextidx = pt2idx(nexttraj.p2) #index of next point (p1 is current pt)
            #print(f'   ... adding {len(self.idxpath)}th pt: {nexttraj.p2} index: {pt2idx(nexttraj.p2)}')

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
        print('2D heuristic path search completed!')
        #print('{:} Total path cost ({:}) = {:8.2f}: '.format(self.searchtype,costtype,self.Tcost))
        #print('idxpath: ', self.idxpath)
        #return path object, float
        self.T_cost = self.cost()  # compute and store traj cost
        #time.sleep(0.5)
        return self, self.T_cost

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



    # plot grid and trajectories (2D only)

    def plot(self,idx, note=''): # plot a path with trajectories
        fig = self.plotSetup(note)
        self.plotOnePath(fig)
        self.plotDone(fig)

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
        arrow_interval = len(self.path) // (Npts-1) # Change 5 to adjust the arrow density
        #print('arrow_interval: {:}  len(self.path) {:}  Npts {:}'.format(arrow_interval, len(self.path), Npts))
        arrow_positions = arrow_positions[:-1:arrow_interval]
        arrow_directions = arrow_directions[::arrow_interval]


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

        # plot the start point as larger dot
        tr = self.path[0]
        #x = [tr.p1.x]
        #y = [tr.p1.v]
        #r = [20]
        startpt = plt.Circle((tr.p1.x,tr.p1.v),0.5,color='green')
        #ax.scatter(x,y,s=r,alpha=0.3)
        ax.add_patch(startpt)

        axlim = 2
        #ax.set_xlim([-axlim,axlim])
        #ax.set_ylim([-axlim,axlim])

    def plotDone(self,figure):
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
            plotSave(figure, my_dpi, imgdir, imgname)

        else:
            print('plot image NOT saved')


class search_from_here6D:
    def __init__(self,Mark,path):
        self.costtype = costtype # a global from config file
        self.pstartIdx = None  # last known trajectory point (or initial point)
        self.cmin   = 99999999999
        self.minidx = 0
        self.minTrs = []
        self.minidxs = []
        self.cminFound = False
        self.ties = None   # number of ties found in this set of branches
        self.mark = Mark  # array to mark already chosen pts (True == still available)
        self.path = path

    def iterate(self,N,function,Cm):
        for ix in range(N):
                for iy in range(N):
                    for iz in range(N):
                        for idx in range(N):
                            for idy in range(N):
                                for idz in range(N):
                                    ivect = [ix,iy,iz,idx,idy,idz]
                                    index = getidx6D(ivect)
                                    function(index,ivect,Cm)
    #6D
    def find_cmin(self,p2idx,ivect,Cm):  # should be called by iterate as first step.
        if self.mark[p2idx]: # index = Cm.m column
            self.IsawOne = True
            p1idx = self.pstartIdx
            if p1idx != p2idx:
                tc = self.eval_cost(p1idx,p2idx,Cm)
                #print('         find_cmin: ', tc, self.cmin)
                ##x = input('  ... pause (CR) ...')
                if tc < self.cmin:
                    self.cmin   = tc
                    self.minidx = p2idx
                    self.cminFound = True

    #6D
    def eval_cost(self,i1,i2,Cm):
        # i1,i2 are indices of traj start and end pt.
        #try:
        #tc,ec = self.path.Cm.m[i1][i2]
        #p1 = pts[i1]
        #p2 = pts[i2]
        #tr = trajectory6D(p1,p2)
        tc,te = getCosts_ij(i1,i2,Cm)
        #except Exception as ex:
            #print('Exception in eval_cost:', type(ex).__name__, ex.args)
            #print('bad path indeces? ',i1,i2)
            #print(f'len pts: {len(pts)}')
            #print(f'examples: 0:{pts[0]}, 9:{pts[9]}')
            ##print('Cm[][]:',self.path.Cm.m[i1][i2])
            #quit()
        if self.costtype == 'energy':
            retCost = te
        elif self.costtype == 'time':
            retCost = tc # precomputed
        else:
            error('search:set_costtype: unknown cost type (6D): '+costtype)
        return retCost

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
        #tr = trajectory6D(point6D(getcoord6D(self.pstartIdx)),point6D(getcoord6D(ti)))
        tr = trajectory6D(pts[self.pstartIdx],pts[ti])
        self.mark[ti] = False  # mark new point as visited
        return ti,tr #chosen next traj index, chosen next traj trajectory

    def find_all_cminTrs(self,index,vect, Cm): # find a list of all next pts for which cost ~= cmin
        epsilon = self.cmin * 0.02 # define 'close'
        if self.mark[index]:
            tc = self.eval_cost(self.pstartIdx,index, Cm)  # get cost for this branch
            if abs(tc-self.cmin) < epsilon:
                tr = trajectory6D(pts[self.pstartIdx],pts[index])
                self.minTrs.append(tr)
                self.minidxs.append(index)
        if len(self.minTrs) > self.path.nmin_max: # find the longest list length
            self.path.nmin_max = len(self.minTrs)
        self.ties = len(self.minidxs) #how many ties at this node??

class path6D:
    #def __init__(self,Cm):
    def __init__(self):
        self.Cm = None     # cost matrix (actually trajectories)
        ##sizer('path6D.Cm: ',self.Cm)
        self.sr = startrow
        self.sc = startcol
        self.mark = [True for x in range(Npts)]  # true if pt is UNvisited
        ##sizer('path6D.mark: ',self.mark)
        self.mark[self.sr*N+self.sc] = False # mark our starting point (can be overridden)
        self.Tcost = 0.0
        self.path = []  # the path as a list of trajectories
        self.idxpath = [] # the path as a list of indices (0..Npts)
        self.searchtype = 'none yet'
        self.datafile = None
        self.maxTiesHSearch = -99999999  # most ties when greedy searching
        self.tie_freq = np.zeros(MAXTIEHISTO)  # histogram of how many ties of each length
        self.ties = 0

    #6D
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
            error('6D grid is to huge for brute force search!')
            #if dfile is None:
                #error('path.search: brute force search requires a dfile')
            #p, cmin = self.bruteForce(dfile=dfile,profiler=profiler)
        #
        # cost-evaluate a random sample of paths
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

    #6D
    def sampleSearch(self,dfile=None,nsamples=977,profiler=None):
        path, cost = self.bruteForce(dfile=dfile,sampling=True,nsamples=nsamples,profiler=profiler)
        return path, cost

    #6D
    def bruteForce(self,dfile=None,sampling=False,nsamples=0,profiler=None): # path class
        if self.searchtype.startswith('none'):
            self.searchtype = 'brute force'
            sampling = False

        if not sampling:
            error('theres NO way I can do full brute force in 6D!!!')

        #n_all_paths = math.factorial(Npts)
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
            tps = [itype]*(Npts)      # path point seq
            tps.append(ftype) # the path cost's type
            names = []
            for i in range(Npts):
                names.append('{:}'.format(i))
            names.append('Cost')
            df.metadata.d['Ncols'] = len(names)
            df.metadata.d['Types'] = tps
            df.metadata.d['Names'] = names
            df.metadata.d['CostType'] = costtype
            df.metadata.d['Space'] = '6D'
            df.metadata.d['Grid dim'] = N
            df.metadata.d['SearchType'] = self.searchtype
            df.metadata.d['#samples'] = nsamples
            #
            df.open()  # let's open the file (default is for writing)

        #
        #  Only sampling - full brute force is deleted
        #
        # 1) generate list of paths
        print('We are generating {:} random sampled paths through {:} nodes'.format(nsamples,Npts))
        phset = set()
        piter = []
        n = 0
        while n < nsamples: # make sure list has no dupes
            if n%10000 == 0:
                print(n,' paths')
            p = list(range(Npts))
            random.shuffle(p) # generate a path as random list of indices
            ph = ''
            for i in range(50):
                ph += str(p[i])[-1] # last digit of idx
            if ph not in phset: # we've found a new pt
                phset.add(ph)
                piter.append(p)
                n+=1 # count adds (faster than len()??)
        #sizer('piter: ',piter)
        print('Path enumeration complete (without duplicates):')

        #2) evaluate their costs

        self.Cm = Cm6D()   # this will store costs for each Tr  as they are encountered
        path_costs = []
        n = -1
        cmin = 99999999999
        cmax = 0
        pmax = path6D()
        pmin = path6D()
        updaterate = 1 + nsamples//20
        for p in piter:  # piter iterates to a series of lists of point indices
            n+=1
            # p is the current path [idx0,idx1,idx2 ...]
            # now get the cost of p
            idxpath = list(p)
            c = cost_idxp6D(costtype, idxpath,self.Cm)  #what is cost of this path?
            if n%updaterate == 0:
                pct = 100.0*n/nsamples
                print(f'seach completion: {pct:4}%') # I'm alive!
            if STOREDATA:
                row = idxpath # list of int index pts
                row.append(c) # tack on the cost
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
                #sizer('pmin: ',pmin)
        #
        #  we are done with the path set to be evaluated
        #
        print(f'{n} paths have been evaluated')
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

    #6D
    def multiHSearch(self,dfile,nsearch,profiler=None):
        ts1 = datetime.datetime.now()
        self.datafile = dfile #keep track of this for adding metadata
        df = dfile
        print('Saving permutations (paths) to: ',df.name)

        # set up the output file metadata
        itype = str(type(5))
        ftype = str(type(3.1415))
        tps = [itype]*(Npts)      # path point-index sequence
        tps.append(ftype) # the path cost's type (float)
        names = []
        for i in range(Npts):
            names.append('{:}'.format(i))
        names.append('Cost')
        df.metadata.d['Types'] = tps
        df.metadata.d['Names'] = names
        df.metadata.d['Ncols'] = len(names)
        df.metadata.d['CostType'] = costtype # 'energy' or 'time'
        df.metadata.d['Space'] = '6D'
        df.metadata.d['Grid dim'] = N
        df.metadata.d['Grid dim'] = N
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
        self.Cm = Cm6D()   # this will store costs for each Tr  as they are encountered

        if nsearch > 2*Npts:  # for big enough searches, allocate same # to all start points
            nperstart = nsearch//Npts
            MULTI_SEARCH_PER_PT = True
            loopcount = Npts
        else:
            nperstart = 1   # if less, just pick random start points
            MULTI_SEARCH_PER_PT = False
            loopcount = nsearch
        print(f'starting {nsearch} 6D iterations with loopcount {loopcount} and nperstart {nperstart}')
        x = input (' ... testpause ...')
        for i in range(loopcount): # go through the start pts
            # go through the start points with equal number at each
            for m in range(nperstart): # do each start pt this many times
                # reset search info
                self.mark = [True for x in range(Npts)]
                count = 0
                if not MULTI_SEARCH_PER_PT:
                    # a random start point
                    startPtIdx = random.randint(0,Npts-1)
                else:
                    startPtIdx = i

                # don't think we need these acth
                self.Tcost = 0.0
                count = 0

                print(f'       6D iteration  {m+1}/{nperstart}  for starting point {i+1}/{loopcount}')
                # do the search
                pself,c = self.heuristicSearch6D(startPtIdx) #including random tie breakers

                if  pself.ties > maxTies:
                    maxTies = pself.ties
                #
                datarow = pself.idxpath
                #
                datarow.append(c) # last col is cost
                df.write(datarow)
                if c < cmin: # find lowest cost of the runs
                    cmin=c
                    pmin = path6D()
                    pmin.path = pself.path
                    pmin.Tcost = c
                if c > cmax: # find highest cost
                    cmax = c
                    pmax = path6D()
                    pmax.path = self.path
                    pmax.Tcost = c

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
            destfolder = df.folder
            if destfolder[-1] != '/':
                destfolder.append('/')
            tfname = destfolder+'ties_info_'+ df.hashcode+'.csv'
            fp = open (tfname, 'w')
            fmtstring1 = '{:8d} , {:8d}, {:8.4f}'
            median=0
            for i,n in enumerate(self.tie_freq):
                median += n
            median /=2
            df.metadata.d['Median Ties']=median
            for i,nst in enumerate(self.tie_freq):
                n = int(nst)
                if i>1:
                    if n == 0:
                        tlv = -1.0 # tmp log value
                    else:
                        tlv = np.log10(n)
                    print(fmtstring1.format(i,n,tlv),file=fp)
            fp.close()
        df.close()
        # return path object, float
        return pmin,pmin.Tcost

    # 6D
    def heuristicSearch6D(self, idx1, profiler=None):
        # sanity check!!
        if Npts > 1.0E4:
            error('too big a search!!: '+float(Npts))

        startPtIdx = idx1  # starting point (<N!)
        self.mark = [True for x in range(Npts)]
        count = 0
        self.mark[idx1] = False # mark our starting point
        self.Tcost = 0.0
        self.path = [pts[startPtIdx]]
        self.idxpath = [startPtIdx] # path has a start point
        self.nmin_max = 0
        latestIdx = startPtIdx
        updaterate = Npts//10
        print(f'                 [..........] path len: {Npts}')
        print(' path completion: ',end='')
        while True:
            li = len(self.path)
            if li%updaterate==0:
                #pct = 1+int(100.0 * li / Npts)
                #print(f' completion: {pct:2}%')
                print('*',flush=True, end='')
            srchFrmHere = search_from_here6D(self.mark,self)
            srchFrmHere.cmin = 99999999 # just being sure/clear
            srchFrmHere.pstartIdx = latestIdx  # updated at end of this loop!
            srchFrmHere.minTrs=[] #these will hold all branches matching cmin cost.
            srchFrmHere.minidxs=[]
            srchFrmHere.IsawOne = False

            ### iteration through the open branches from this node
            srchFrmHere.iterate(N,srchFrmHere.find_cmin, self.Cm)
            #DEBUGgate = li > 725 or li < 4
            DEBUGgate = False
            if DEBUGgate:
                print(f'pL={li}: Started at: {latestIdx}: any mins?: {srchFrmHere.IsawOne}, minfound: {srchFrmHere.cminFound}, mincost: {srchFrmHere.cmin}')
            if not srchFrmHere.cminFound:
                error('heuristic search iteration: somethings wrong, need to find_cmin before find_all')
            srchFrmHere.iterate(N,srchFrmHere.find_all_cminTrs, self.Cm)
            ###

            # keep track of greatest number of ties along this path
            if srchFrmHere.ties > self.maxTiesHSearch:
                self.maxTiesHSearch = srchFrmHere.ties

            if DEBUGgate:
                print(f'pL={li}:     Ties: {len(srchFrmHere.minidxs)}/{len(srchFrmHere.minTrs)}')
            nxtidx,nxtTr = srchFrmHere.select_next()   # break a possible tie btwn branches leaving this pt.
            if DEBUGgate:
                print(f'pL={li}:  ...   appending {nxtidx}')
                if nxtidx in self.idxpath:
                    print(f'      {idxpath} is already in!')
            self.path.append(nxtTr)
            if DEBUGgate:
                print(f'pL={li}:      after append, len(self.path)={len(self.path)}')
            self.idxpath.append(nxtidx)
            latestIdx = nxtidx
            if DEBUGgate:
                print(f'pL={li}:       exit condition: {len(self.path)},{Npts}-> {len(self.path) >= Npts}')
            if len(self.path) >= Npts:  # pts path has the full Npts (not traj's)
                break
        print('')
        self.Tcost = cost_idxp6D(costtype, self.idxpath,self.Cm)  # compute total cost of the path
        print('heuristic path search completed!')
        print('{:} Total path cost ({:}) = {:8.2f}: '.format(self.searchtype,costtype,self.Tcost))
        # return path object, float
        return self, self.Tcost
    # 6D
    def check(self): # 6D
        if len(self.path) != Npts-1:
            error('wrong path length '+str(len(self.path)))
        i=0
        for t in self.path:
            if not t.valid:
                error('path contains invalid traj: '+str(i))
            if i > 0:
                if self.path[i-1].p2 != t.p1:
                    error('discontinuous path: '+str(i-1)+' '+str(i))
            i+=1
    #6D
    def compute_curves6D(self,idx): #6D
        curvepts_x = []
        for i,tr in enumerate(self.path):
            if i in idx:
                if i == len(self.path):
                    break
                if not tr.valid:
                    error('compute_curves()6D: I should not have found an invalid trajectory in path: '+str(i))
                if tr is None:
                    error('null traj: '+ str(i) + str(tr))
                if not tr.computed and not tr.constrained:
                    error('Cant plot until trajectory is computed and constrained '+ str(i) + str(tr))
                dt = tr.dt
                for i in range(NPC):
                    t = dt*i/NPC
                    curvepts_x.append(tr.x(t))
                curvepts_x.append(tr.x(dt))
        x = np.array(curvepts_x).T
        return x

    #6D
    def save6D(self,hashcode): # save 6D trajectory for animation and plotting
        #
        #   this is a new datafile just for visualization
        #
        df = bd.datafile('6Dtrajdata', 'BH', 'simulation')
        df.hashcode = hashcode # keep hashcode same as search df.
        df.set_folders('','') # default local folders
        trajcurves = self.compute_curves6D(-1) # save all trajectories
        print('path.save: x points:     ',len(self.path))
        print('path.save: x curves dims:', trajcurves.shape)
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
        r,c = np.shape(trajcurves)
        for i in range(c):
            row = [i,tracurves[0][i],tracurves[1][i],tracurves[2][i]]
            df.write(row)

        df.close()   # all done


def main():
    import pickle
    import os

    epsilon = 0.02

    configure()

    print('Commencing tests: 6D, N=',N)

    print('index <--> coordinates test')

    v=[0,0,0,0,0]  # 5 random 6D test coordinate vectors

    for i in range(len(v)):
        v[i] = [0]*6
        for j in range(6):
            v[i][j] = random.choice(range(N))
        print('v: [',i,']:',v[i])

    for i in range(3):
            x = getcoord6D(getidx6D(v[i]))
            assert x==v[i]

    print('index <--> coordinates test: PASSED')

    # get curves for each arc
    cx, cy = self.compute_curves(-1) #compute trajectory path
    ax.plot(cx,cy,color='blue')
    axlim = 2
    ax.set_xlim([-axlim,axlim])
    ax.set_ylim([-axlim,axlim])

    print('\n\n   Cm tests: ')

    c1 = Cm2D()

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
    assert r == Npts

    p1 = point6D(getcoord6D(1234))
    p2 = point6D(getcoord6D(3241))
    tr12 = trajectory6D(p1,p2)
    tt3d = type(tr12)

    if not SKIPFILL:
        # [10][10] is just a "random" traj
        assert type(c1.m[10][10]) == tt3d
        assert not c1.m[10][10].valid  # self-self transitions invalid.
        assert c1.m[10][11].valid
    else:
        print('\n                 fill test skipped!!\n')
        c1.m[10][10] = tr12

    # try a path:
    p = path2D(gt,c1,r,c)
    p.heuristicSearch()

    print(p)
    print('   Cm tests:  PASSED')

    print('\n\n   Amax tests')
    print('\n\n If this test fails, try a smaller value of DT_START (in ctoConfig.txt)\n')

    print('DT_START: ',DT_START)
    print('AMAX: ',AMAX)
    for i in range(10):
        p1 = point6D(getcoord6D(random.randint(0,Npts-1)))
        p2 = point6D(getcoord6D(random.randint(0,Npts-1)))
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
    tr12.getCosts()
    #a = tr12.timeEvolution(ACC_ONLY=True)
    #tr12.cost_e(a)
    #tr12.cost_t(a)
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
    bpath = path6D()
    bpath.search(searchtype,dfile=df,nsamples=10000)


    print('\n\n            ALL tests:     PASSED')

if __name__ ==  '__main__':
    print('main starting:')
    main()

