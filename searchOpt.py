#!/usr/bin/python3

import logging
import c2to as cto
import sys
#import matplotlib.pyplot as plt


def main(args):
    SEARCHT = 'heuristic search'
    #SEARCHT = 'brute force'

    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    # create grid
    gt = cto.grid2D(cto.N)

    # cost matrix
    c1 = cto.Cm()
    c1.fill(gt)


    # compute a path:
    p = cto.path(gt,c1)

    if SEARCHT.startswith('heur'):
        path2, cmin = p.heuristicSearch()
    elif SEARCHT.startswith('brute'):
        path2, cmin = p.bruteForce(c1)

    # is it a valid path?
    #p.check()

    notes = '{:}, cost: {:4.2f} ({:})'.format(SEARCHT, cmin, cto.costtype)

    # graph the path
    path2.plot(idx,notes)
    #p.plot(idx)





if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
