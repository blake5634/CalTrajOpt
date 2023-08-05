#!/usr/bin/python3

import logging
import c2to as cto
# these are needed for unpickling
#from c2to import Cm as Cm
#from c2to import trajectory3D as trajectory3D
#from c2to import point3D as point3D
import sys
import os
import pickle
import brl_data.brl_data as bd
import datetime as dt
import random

#import matplotlib.pyplot as plt


def main(args):

    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    #c1 = cto.Cm(df = df)  # cost matrix



    ##########################################################
    #
    #  Grid Type:
    RANDOMGRID = True
    N = 4
    cto.N = N

    # cost:
    cts = ['time','energy']
    #cto.costtype = 'energy'
    cto.costtype = 'time'
    #
    ###########################################################
    if RANDOMGRID:
        #df.metadata.d['Random Grid'] = True
        cto.gridType = 'random'
    else:
        #df.metadata.d['Random Grid'] = False
        cto.gridType = 'rectangular'

    pl = list(range(N**6))
    random.shuffle(pl) # generate a path as random list of indices
    for i in range(10):
        print(i+500, pl[i+500])
    #
    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D()
    p.idxpath = pl

    pts = cto.setupPoints()   # just store points instead of cost matrix Cm
    L = len(p.idxpath)-1
    print('\n\n')
    print('length of point list: ',len(pts))

    if L != N**6-1:
        error('not a correct length path!: '+str(L))

    assert L+1 == len(pl),  'L problem'
    print('max(pl): ',max(pl))
    assert len(pl)-1 == max(pl), 'max problem'
    print('path len:',len(p.idxpath))
    connectedCt = 0
    pprev = pts[p.idxpath[0]]
    for i in range(L):
        i1 = p.idxpath[i]
        i2 = p.idxpath[i+1]
        p1 = pts[i1]
        p2 = pts[i2]
        cs = 'no'
        if p1.xvect == pprev.xvect:
            connectedCt +=1
            cs = 'yes'
        if i1==i2:
            error('path repeats a node index')
        #print( i,i+1)
        #print('\n ', i1,i2)
        #print(pts[i1],'  -->  ',pts[i2], cs)
        pprev = p2
        #x = input('pause ...')
        #ct,ce = Cm.m[i1][i2]
    #c = cto.cost_idxp(costtype, p.idxpath)  #what is cost of this path?
    assert connectedCt == len(p.idxpath)-1

    print('\n\n                      3D/6D path generation continuity PASSED \n\n')

if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
