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

def approx(x,y,epsilon):
    return abs((x-y)/y) < epsilon

if DIMENSION not in ['2D','6D']:
    cto.error("DIMENSION must be in ['2D', '6D']")
else:
#def main(args):
    # create a path and plot it graphically
    cto.configure()
    ###################################################
    #    Search type and size  (some of these matter to points storage testing)
    #

    gridtype = 'random'
    #gridtype = 'rectangular'

    cto.N = 4
    cto.M = 16
    N = cto.N
    M = cto.M
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


    # create a datafile to store  random points

    df = bd.datafile('randomGridPointSet','BH','simulation')
    dataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom'
    codeFolder = ''
    df.set_folders(dataFolder,codeFolder) # also creates filename


    #
    # create grid and cost matrix
    if DIMENSION == '2D':
        gt = cto.grid2D(cto.N)
        c1 = cto.Cm()
        if gridtype == 'random':
            gt.randgrid = True
            #  sadly, grid is randomize through the Cost matrix (Cm) fill() method!
            c1.set_GridRandomize(df=df)  # select random instead of grid

        c1.fill(gt) # optionally randomize grid and fill up costs in Cm.m[][]
        x = input('    pause ...')
    else:
        cto.error('6D not yet implemented')


    # generate the points
    if DIMENSION == '2D':
        gt.savePoints2D(df) #  NEW ( a grid method )

        # store the points which were ostensibly saved:
        grid_saved = cto.grid2D(cto.N) # new 2D grid
        for i1 in range(N):
            for j1 in range(N):
                grid_saved.gr[i1][j1] = gt.gr[i1][j1] # not cost matrix!

        print(f"saved random point set to {df.name}.")


        # read the points back
        gt.readPoints2D(df)  # NEW ( a grid method)
        print(f'points read in from {df.name}')

        # compare coordinates of saved vs. read points
        for i1 in range(N):
            for j1 in range(N):
                xmatch = approx(grid_saved.gr[i1][j1].x,gt.gr[i1][j1].x,0.0002)
                vmatch = approx(grid_saved.gr[i1][j1].v,gt.gr[i1][j1].v,0.0002)
                assert xmatch and vmatch , f'storage mismatch: [{i1}][{j1}]'

        print('\n\n         2D Point store/read operation test          PASSED')
        quit()

    else: #6D
        cto.error('6D not yet implemented')

