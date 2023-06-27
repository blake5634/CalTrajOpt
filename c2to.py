#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import brl_data.brl_data as bd
import sys
import random

def error(msg):
    print('Error: ')
    print(msg)
    quit()


global costtype, AMAX, DT_TEST, N,M, NPC, startrow, startcol

N = 4
M = N*N  # number of distances

AMAX = 2  #  normalize for now.
DT_TEST = 2.0


NPC = 20  # number plot points per curve

costtype = 'time'
costtype = 'energy'

# for 1D optimization case
startrow = 3
startcol = 1

def configure(fp=None):
    #print('Pyton path: ', sys.path)
    global costtype, AMAX, DT_TEST, N,M, NPC, startrow, startcol
    if not fp:
        f = open('ctoConfig.txt','r')
    else:
        f = fp

    for line in f:
        parname, *rest = line.split(' ')
        parname = parname.strip()
        v = ' '.join(rest).strip() # parameter value as string

        if parname == '#':
            next

        if parname =='costtype':
            if v not in ['time','energy']:
                error('path.configure: unknown cost type: '+v)
            costtype = v
        if parname == 'amax':
            AMAX = float(v)
        if parname == 'dt_test':
            DT_TEST = dt_test  # fixed dt value used for testing
        if parname == 'N':
            N = int(v)
            M = N*N  # used in 1D code
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

    def __init__(self,iv):
        if len(iv) != 6:
            error('point3D: invalid index vector')
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
        if str(type(p1)) != "<class 'c2to.point3D'>" or str(type(p2)) != "<class 'c2to.point3D'>":
            error('trajectory3D called with 1D point!')
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
        dt = 1.4
        ni = 0
        while True:
            ni += 1
            self.compute(dt)
            if self.get_Amax(dt) > AMAX:
                dt *= 1.05
            else:
                break
        dt *= 0.95
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
        return c

    def cost_t(self, a):  # need a parameter for efficient calling though not used(!)
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

    def timeEvolution(self,ACC_ONLY = False):
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
            if not ACC_ONLY:
                x.append(self.x(t2))
                v.append(self.xd(t2))
            a.append(self.xdd(t2))
        if not ACC_ONLY:
            t.append(self.dt)
            x.append(self.x(self.dt))
            v.append(self.xd(self.dt))
        a.append(self.xdd(self.dt))
        if ACC_ONLY:
            return a
        else:
            return t,x,v,a

    #
    #   Energy cost: sum of squared acceleration
    def cost_e(self,a):
        if not self.computed and not self.constrained:
            error('Cant compute cost_e until trajectory is constrained')
        c = 0.0
        n = len(a)
        for a1 in a:
            c += self.dt*a1*a1/n
        return c

    def cost_t(self):
        if not self.computed and not self.constrained:
            error('Cant compute cost_t until trajectory is constrained')
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

class search_from_curr_pt:
    def __init__(self,Mark):
        self.costtype = 'None yet'
        self.pstart = None  # last known trajectory point (or initial point)
        self.cmin   = 99999999999
        self.minidx = 0
        self.minTrs = []
        self.minidxs = []
        self.found = False
        self.mark = Mark  # array to mark already chosen pts (True == still available)


    def iterate(self,N,function):
        def getidx(v):
            return v[0]*N**5+v[1]*N**4+v[2]*N*N*N+v[3]*N*N+v[4]*N+v[5]
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
        if self.mark[index]:
            p1 = self.pstart
            p2 = point3D(ivect)
            if p1==p2:
                p2.valid = False
            if p2.valid:
                tr = trajectory3D(p1,p2)
                tc = self.eval_cost(tr)
                if tc < self.cmin:
                    self.cmin   = tc
                    self.p2Cmin = tr.p2
                    self.trCmin = tr
                    self.minidx = index
                    self.found = True

    def eval_cost(self,tr):
        tr.constrain_A()
        if self.costtype == 'energy':
            a =  tr.timeEvolution(ACC_ONLY=True)
            tc = tr.cost_e(a)
        elif self.costtype == 'time':
            tc = tr.cost_t()
        else:
            error('search:set_costtype: unknown cost type (3D): '+costtype)
        return tc

    def select_next(self):
        L = len(self.minTrs)
        if L < 1:
            error('search.select_next: no next trajs identified yet. did you run find_all_cminTrs?')
        if L == 1:
            self.mark[self.minidx] = False
            return self.minidx, self.minTrs[0]
        # pick one at random
        print('select_next: choosing a random min traj! from ',len(self.minidxs))
        tr = random.choice(self.minTrs)   #
        if not tr.computed or not tr.constrained:
            error('select_next: trajectory is not comp or const.')
        ti = random.choice(self.minidxs)
        self.mark[ti] = False
        return ti,tr

    def find_all_cminTrs(self,index,ivect): # find a list of all next pts for which cost ~= cmin
        if not self.found:
            error('search.find_all_cminTrs: somethings wrong, need to find_cmin before find_all')
        epsilon = self.cmin * 0.05 # define 'close'
        if self.mark[index]:
            p1 = self.pstart
            p2 = point3D(ivect)
            if p1==p2:
                p2.valid = False
            if p2.valid:
                tr = trajectory3D(p1,p2)
                tc = self.eval_cost(tr)
                if abs(tc-self.cmin) < epsilon:
                    self.minTrs.append(tr)
                    self.minidxs.append(index)


# find the path without filling up the cost matrix Cm
class path3D:
    def __init__(self, adv=True):
        ADVANCED = adv
        BASIC = not ADVANCED
        self.advflag = adv

        # sanity check!!
        if N**6 > 1.0E4:
            error('too big a search!!: '+float(N**6))
        def getidx(v):
            return v[0]*N**5+v[1]*N**4+v[2]*N*N*N+v[3]*N*N+v[4]*N+v[5]

        v1 = [1,2,1,1,2,0]
        pstart = point3D(v1)  # starting point (<N!)
        mark = [True for x in range(N**6)]
        count = 0
        mark[getidx(v1)] = False # mark our starting point

        if ADVANCED:
            self.Tcost = 0.0
            self.path = []
            plen = 0
            nmin_max = 0
            while len(self.path) < N**6 - 1:
                search = search_from_curr_pt(mark)
                search.costtype = costtype
                search.pstart = pstart
                search.iterate(N,search.find_cmin)
                search.iterate(N,search.find_all_cminTrs)
                if len(search.minidxs) > nmin_max:
                    nmin_max = len(search.minidxs)
                ni,nr = search.select_next()
                self.path.append(nr)
                plen += 1
                pstart = nr.p2
                self.Tcost += search.cmin

            print('\n\n       The longest set of min-cost next points was: {:} points\n\n'.format(nmin_max))
            print('ADVANCED Path search completed!')
            print('Total path cost ({:}) = {:8.2f}: '.format(costtype,self.Tcost))

        if BASIC:
            self.Tcost = 0.0
            self.path = []
            plen = 0
            while len(self.path) < N**6-1:
                # start adding next point to path list
                cmin = 99999999
                found = False
                # search costs through the available next points:
                for ix in range(N):
                    for iy in range(N):
                        for iz in range(N):
                            for idx in range(N):
                                for idy in range(N):
                                    for idz in range(N):
                                        index = getidx([ix,iy,iz,idx,idy,idz])
                                        if mark[index]:
                                            p1 = pstart
                                            p2 = point3D([ix,iy,iz,idx,idy,idz])
                                            if p1==p2:
                                                p2.valid = False
                                            if p2.valid:
                                                #
                                                #  simple greedy search:
                                                #       take the first min value encountered
                                                #       (multiple points may have ~= cost
                                                #
                                                tr = trajectory3D(p1,p2)
                                                tr.constrain_A()
                                                if costtype == 'energy':
                                                    a =  tr.timeEvolution(ACC_ONLY=True)
                                                    tc = tr.cost_e(a)
                                                elif costtype == 'time':
                                                    tc = tr.cost_t()
                                                else:
                                                    error('unknown cost type (3D): '+costtype)
                                                if tc < cmin:
                                                    cmin = tc
                                                    p2Cmin = p2
                                                    trCmin = tr
                                                    minidx = index
                                                    found = True

                # we've found lowest cost cmin and that point
                #
                # make sure though:
                if not found:
                    print('current index: ',index,'/', N**6)
                    error('I didnt find a minimum cost!!'+ str(count))
                if not mark[minidx]:
                    error('new path point is marked already!')
                if not p2.valid:
                    print('current index: ',index,'/', N**6)
                    error('path: trying to add invalid new point' + str(p2))
                # future: go back and find all points with same cost as trCmin
                if plen > 0:
                    pprev = self.path[-1].p2
                    pcurr = trCmin.p1
                    if pprev != pcurr:
                        print('minidx: ', minidx)
                        print('adding traj: ', trCmin)
                        error('path trajectories dont connect! '+str(pprev))
                # we've passed all the checks: append the traj to path
                print('adding: ', plen, trCmin.p2)
                mark[minidx] = False  # do not visit this point again
                self.path.append(trCmin)
                plen += 1
                pstart = p2Cmin
                self.Tcost += cmin
            print('BASIC Path search completed!')
            print('Total path cost ({:}) = {:8.2f}: '.format(costtype,self.Tcost))

    def check(self): # 3D
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


    def save(self,fname): # 3D
        df = bd.datafile(fname, 'BH', 'simulation')
        df.set_folders('','') # default local foldesr
        x = self.compute_curves(-1) # save all trajectories
        print('path.save: x dims:', x.shape)
        col_names = ['n','X','Y','Z']
        int_type = str(type(5))            # these have to be strings b/c json can't serialize types(!)
        float_type = str(type(3.14159))
        col_types = [ int_type, float_type, float_type, float_type]
        col_comments = ['','','','']
        df.set_metadata(col_names, col_types, col_comments)
        df.metadata.d['AMAX']=AMAX
        df.metadata.d['N'] = N
        df.metadata.d['costtype'] = costtype
        df.metadata.d['Advanced_searchtype'] = self.advflag
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
                        a = Cm.m[crow][ccol].timeEvolution(ACC_ONLY=True)
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

