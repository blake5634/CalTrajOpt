#!/usr/bin/python3

import c2to as cto
import sys
import os
import pickle
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt

def fill_check(Cmat):
    print('checking the fill...')
    nf = 0
    M = cto.M
    for i1 in range(M):  # go through all grid points
        for j1 in range(M):
            if Cmat.m[i1][j1] == 0:
                print('{:},{:} has zero value'.format(i1,j1))
            if type(Cmat.m[i1][j1]) != type((5.1,6.2)):
                print('{:},{:}: has type {:}'.format(i1,j1,type((5.1,6.2))))
            if i1 ==0 and j1 == 364:
                ct,ce = Cmat.m[i1][j1]
                print('{:},{:}: ct:{:4.1f} ce:{:4.1f}'.format(i1,j1,ct,ce))
            if i1 ==364 and j1 == 364:
                try:
                    ct,ce = Cmat.m[i1][j1]
                    print(' ... i1 == j1  worked')
                except Exception as extype:
                    print('Error: {:}, at {:},{:}: ct:{:4.1f} ce:{:4.1f}'.format(extype.__name__,i1,j1,-999.9, -999.9))
    print('checking the fill COMPLETED')


def main(args):
    q = 'debug Cm filling'

    # create the datafile:
    df = bd.datafile('6Dsearching','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')

    df.metadata.d['Research Question'] = q

    # create a path and plot it graphically
    cto.configure()

    c1 = cto.Cm(df = df)  # cost matrix
    c1.set_GridRandomize()  # needed? redundant??
    c1.fill()

    fill_check(c1)

    if False:

        ##########################################################
        #
        #  Grid Type:
        RANDOMGRID = True   # remove the pickle file when changed!!!

        #
        # 1D only
        #searchtype = 'brute force'
        # 4x4 and 6D:
        searchType = 'multi heuristic'
        #searchType = 'sampling search'

        # approx 1M samples
        nsamp = 4*cto.N**6  # at least one for each starting point(!)
        #nsamp = N**6 # 1M

        #nsamp = 60 # for testing

        # cost:
        cts = ['time','energy']
        #cto.costtype = 'energy'
        cto.costtype = 'time'
        #
        ###########################################################

        #if RANDOMGRID:
            #df.metadata.d['Random Grid'] = True
        #else:
            #df.metadata.d['Random Grid'] = False

        # memory profiling
        mem_snap('populate Cm')

        pname = 'c1Costs.pickle'
        if os.path.exists(pname):
            print('loading precomputed cost matrix   ...')
            f = open(pname, 'rb')
            c1 = pickle.load(f)
            if RANDOMGRID:
                c1.set_GridRandomize()  # needed? redundant??
            print(' loading completed')
        else:
            print('no stored data: computing cost matrix')
            if RANDOMGRID:
                #print('    ... bug: always have to recompute in Random Grid mode!!!')
                c1.set_GridRandomize()  # random pts instead of grid
            c1.fill()
            f = open(pname,'wb')
            pickle.dump(c1,f)


    print('Completed:')

    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    #p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
