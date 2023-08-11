#!/usr/bin/python3

import numpy as np
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt


def main(args):
    cto.configure()
    ###################################################
    #    Search grid, type and size
    #

    # this code can do two things: generate grid, search grid.
    POINTS_MODE = 'search'
    #POINTS_MODE = 'generate'
    gridtype = 'random'
    ##gridtype = 'rectangular'
    DataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    #pointsDataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom'
    #pointsDataName = '2023-08-11_0f616e46_randomGridPointSet_BH_simulation.csv'  # for reading

    #pointsFilename = pointsDataFolder + '/'+pointsDataName  # full path of points file

    myf = bd.finder()
    myf.set_dirs(['/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom'])
    keys = [args[1], '.csv']
    fns = myf.findh(keys)
    pointsDataFolder = fns[0][0]
    pointsDataName   = fns[0][1]

    pointsFilename = pointsDataFolder + '/'+pointsDataName  # full path of points file

    cto.N = 3
    cto.M = cto.N*cto.N
    #SEARCHT = 'heuristic search' # greedy nearest neighbor
    #SEARCHT = 'brute force'   # enumerate all paths
    #SEARCHT = 'sampling search' # nsearch random paths
    SEARCHT = 'multi heuristic' # repeated heuristic search all starting pts

    #nsearch = int(np.math.factorial(9)* 0.10)  # 10% of 3x3
    #nsearch = 1000000  # 1M
    #nsearch = 4*cto.M   # 4 searches from each starting pt
    nsearch = 10
    #cto.costtype = 'time'
    cto.costtype = 'energy'
    cto.NPC = 30   #  # of simulation points in 0-dt time intervale
    #
    ###################################################
    if POINTS_MODE not in ['generate', 'search']:
        cto.error('searchOpt.py: illegal POINTS MODE')

    if POINTS_MODE == 'generate' and gridtype != 'random':
        cto.error('SearchOpt: Cannot generate points file unless grid is random')

    #idx = int(args[1])
    idx = -1

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    c1 = cto.Cm()
    if gridtype == 'random':
        gt.randgrid = True
        c1.set_GridRandomize()  # select random instead of grid

    if POINTS_MODE == 'generate':
        df = bd.datafile('randomGridPointSet','BH','simulation')
        dataFolder = pointsDataFolder
        codeFolder = ''
        df.set_folders(dataFolder,codeFolder) # also creates filename
        c1.fill(gt)  # randomize
        gt.savePoints2D(df)
        print(f'Random pts saved to {df.name}')
        notes = f'generated random points file: {cto.N}x{cto.N}'
        logentry(df,notes)
        print('\n\n               your points data file hash is:',df.hashcode)

    else: # search mode with grid or read-in points
        if gridtype == 'random': # we should read from file for repeatability
            # create the random points file reader
            dfr = bd.datafile('','','') #'' ok for reading
            dfr.set_folders('','') # '' ok for reading
            dfr.name = pointsFilename
            pointSourceFileName = gt.readPoints2D(dfr)  #read in the set of random points
        c1.fill(gt) # calc trajectories and costs after points reading
        dfw = bd.datafile('2Dsearching','BH','simulation')
        dfw.set_folders(DataFolder,'')
        dfw.metadata.d['Points Data Source'] = pointSourceFileName
        q = input('Research Question for this search:')
        dfw.metadata.d['Research Question'] = q
        # instantiate a path:
        p = cto.path(gt,c1)
        path2, cmin = p.search(SEARCHT, dfile=dfw, nsamples=nsearch)
        print('Optimal path returned: (tra)', path2.path)
        print('Optimal path returned: (idx)', path2.idxpath)
        # is it a valid path?
        #p.check()
        notes = f"{gridtype} grid, {SEARCHT}, cost: {cmin:8.1f} ({cto.costtype})"   #  keep a "log book"
        logentry(dfw,notes)
    #graph the path
    if POINTS_MODE == 'search':
        path2.plot(-1,notes)
    #p.plot(idx)

def logentry(df,notes):
    logdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    logfilename = 'work_logbook.txt'
    q = df.metadata.d['Research Question']
    if len(q)>0 and 'debug' not in q:  # skip junk files
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        logentry = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, q)
        f =open(logdir+logfilename,'a')
        print(logentry,file=f)
        f.close()
    else:
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')

if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
