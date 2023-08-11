#!/usr/bin/python3

import numpy as np
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt


def main(args):
    #df = None  # no data file output
    df = bd.datafile('path_search_results','BH','simulation')
    dataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    codeFolder = ''
    df.set_folders(dataFolder,codeFolder) # also creates filename

    if df is not None:
        df.metadata.d['ResearchQuestion'] = input('Enter research question:')
    #  cto.point2D.search() will take care of metatada setup

    # create a path and plot it graphically
    cto.configure()
    ###################################################
    #    Search type and size
    #
    gridtype = 'random'
    #gridtype = 'rectangular'

    cto.N = 3
    cto.M = cto.N*cto.N
    #SEARCHT = 'heuristic search' # greedy nearest neighbor
    #SEARCHT = 'brute force'   # enumerate all paths
    #SEARCHT = 'sampling search' # nsearch random paths
    SEARCHT = 'multi heuristic' # repeated heuristic search all starting pts

    #nsearch = int(np.math.factorial(9)* 0.10)  # 10% of 3x3
    #nsearch = 1000000  # 1M
    nsearch = 4*np.math.factorial(cto.M)  # 4 searches from each starting pt

    #cto.costtype = 'time'
    cto.costtype = 'energy'
    cto.NPC = 30   #  # of simulation points in 0-dt time intervale
    #
    ###################################################

    #idx = int(args[1])
    idx = -1

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    c1 = cto.Cm()
    if gridtype == 'random':
        c1.set_GridRandomize(df=df)  # select random instead of grid

    c1.fill(gt)
    x = input('    pause ...')


    # instatntiate a path:
    p = cto.path(gt,c1)
    path2, cmin = p.search(SEARCHT, dfile=df, nsamples=nsearch)
    print('Optimal path returned: (tra)', path2.path)
    print('Optimal path returned: (idx)', path2.idxpath)
    # is it a valid path?
    #p.check()

    notes = f"{gridtype} grid, {SEARCHT}, cost: {cmin:8.1f} ({cto.costtype})\n      ({df.hashcode})"

    #  keep a "log book"

    logdir = dataFolder+'/writing/'
    logfilename = 'work_logbook.txt'

    q = df.metadata.d['ResearchQuestion']
    if len(q)>0 and 'debug' not in q:  # skip junk files
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        logentry = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, q)
        f =open(logdir+logfilename,'a')
        print(logentry,file=f)
        f.close()
    else:
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')

    print('\n\n               your data file hash is:',df.hashcode)
    #graph the path
    path2.plot(-1,notes)
    #p.plot(idx)


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
