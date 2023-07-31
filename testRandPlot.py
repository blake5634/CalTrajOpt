#!/usr/bin/python3

import numpy as np
import math
import matplotlib.pyplot as plt
import itertools as itt
import random
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt

def ru():
    return random.uniform(-1,1)

def main(args):

    RANDOM = True

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    c1 = cto.Cm()
    c1.set_GridRandomize()
    c1.fill(gt)

    # instantiate a path:
    p = cto.path(gt,c1)
    p.idxpath = []
    tv = []
    for i in range(16):
        tv.append(i)
    print(tv)
    random.shuffle(tv) # inplace!!
    p.idxpath = list(tv)
    print(p.idxpath)

    for i in range(len(p.idxpath)-1):
        # build next trajectory
        row,col = cto.idx2ij(p.idxpath[i])
        p1 = cto.point2D(row,col)
        row,col = cto.idx2ij(p.idxpath[i+1])
        p2 = cto.point2D(row,col)
        tr = cto.trajectory2D(p1,p2)
        if RANDOM:
            tr.p1.x = ru()
            tr.p1.v = ru()
            tr.p2.x = ru()
            tr.p2.v = ru()
        tr.constrain_A()
        p.path.append(tr)

    p.plot(8,'testing testing')


    ptidx = 8  # th traj in the path
    r,c = cto.idx2ij(p.idxpath[ptidx+1])
    tt = p.Cm.m[r][c]
    tt.constrain_A()
    print('self.idxpath stored path 8traj:')
    t1 = p.path[ptidx]
    print('start: {:4.2f}, {:4.2f}'.format(t1.p1.x,t1.p1.v))
    print('goal:  {:4.2f}, {:4.2f}'.format(t1.p2.x,t1.p2.v))
    print('Cm stored path 8 traj:')
    print('start: {:4.2f}, {:4.2f}'.format(tt.p1.x,tt.p1.v))
    print('goal:  {:4.2f}, {:4.2f}'.format(tt.p2.x,tt.p2.v))
    print('traj: ',tt)
    print('start/goal x,xd: ({:4.2f},{:4.2f}) --> ({:4.2f}, {:4.2f})'.format(tt.x(0),tt.xd(0),tt.x(tt.dt),tt.xd(tt.dt)))
    p.plotDone()


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
