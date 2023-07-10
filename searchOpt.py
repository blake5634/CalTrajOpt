#!/usr/bin/python3

import logging
import c2to as cto
import sys
import brl_data.brl_data as bd
#import matplotlib.pyplot as plt


def main(args):
    #SEARCHT = 'heuristic search' # greedy nearest neighbor
    SEARCHT = 'brute force'   # enumerate all paths
    #SEARCHT = 'sampling search' # random paths
    SEARCHT = 'multi heuristic' # repeated heuristic search all starting pts

    # create a path and plot it graphically
    cto.configure()
    #idx = int(args[1])
    idx = -1

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    c1 = cto.Cm()
    c1.fill(gt)


    # compute a path:
    p = cto.path(gt,c1)

    #df = None  # no data file output
    df = bd.datafile('path_search_results','BH','simulation')
    dataFolder = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    codeFolder = ''
    df.set_folders(dataFolder,codeFolder) # also creates filename

    if df is not None:
        df.metadata.d['ResearchQuestion'] = input('Enter research question:')
    #  cto.point2D.search() will take care of metatada setup


    ###################################################
    #    Search Size
    #
    nsearch = 100000
    #
    ###################################################

    path2, cmin = p.search(SEARCHT, dfile=df, nsamples=nsearch)

    # is it a valid path?
    #p.check()

    notes = '{:}, cost: {:4.2f} ({:})'.format(SEARCHT, cmin, cto.costtype)

    print('your data file hash is:',df.hashcode)
    # graph the path
    #path2.plot(idx,notes)
    #p.plot(idx)


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
