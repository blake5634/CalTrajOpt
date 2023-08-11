#!/usr/bin/python3

import numpy as np
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt

#####################################################################
#
#   Test saving and reading (randomized) points to/from a file
#      purpose: allow multiple search runs on the SAME random points
#
#     2D case: branch main
#     6D case; branch multiOpt
#          (this code should work on both braches if line 21 is set
#           correctly)
####################################################################

DIMENSION = '2D'  # or '6D'


if DIMENSION not in ['2D','6D']:
    cto.error('DIMENSION must be in ['2D', '6D']')

#def main(args):
    # create a path and plot it graphically
    cto.configure()
    ###################################################
    #    Search type and size  (some of these matter to points storage testing)
    #

    #gridtype = 'random'
    gridtype = 'rectangular'

    cto.N = 4
    cto.M = 16
    #SEARCHT = 'heuristic search' # greedy nearest neighbor
    #SEARCHT = 'brute force'   # enumerate all paths
    #SEARCHT = 'sampling search' # nsearch random paths
    #SEARCHT = 'multi heuristic' # repeated heuristic search all starting pts

    #nsearch = int(np.math.factorial(9)* 0.10)  # 10% of 3x3
    #nsearch = 1000000  # 1M
    #nsearch = 5

    #cto.NPC = 30   #  # of simulation points in 0-dt time intervale
    #
    ###################################################

    #idx = int(args[1])
    idx = -1

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    if DIMENSION == '2D':
        c1 = cto.Cm()
        if gridtype == 'random':
            c1.set_GridRandomize(df=df)  # select random instead of grid

        c1.fill(gt)
        x = input('    pause ...')
    else:
        cto.error('6D not yet implemented')

    # create a datafile to store teh random points

    df = bd.datafile('randomGridPointSet','BH','simulation')
    dataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    codeFolder = ''
    df.set_folders(dataFolder,codeFolder) # also creates filename
    df.metadata.d['ResearchQuestion'] = 'RandomPointSet'

    # generate the points
    if DIMENSION == '2D':
        c1.savePoints2D(df)
        # store the points which were ostensibly saved:
        cm_save = cto.Cm()
        for i1 in range(N):
            for j1 in range(N):
                cm_save.m[i1][j1] = c1.m[i1][j1]

        # read the points back
        c1.readPoints2D(df)

        # compare coordinates of saved vs. read points
        for i1 in range(N):
            for j1 in range(N):
                assert cm_save.m[i1][j1] == c1.m[i1][j1], f'storage mismatch: [{i1}][{j1}]'


    else:
        cto.error('6D not yet implemented')

    print('\n\n         Point store/read operation test          PASSED')
