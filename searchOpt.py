#!/usr/bin/python3

import logging
import c2to as cto
import sys
#import matplotlib.pyplot as plt


def main(args):
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
    p.heuristicSearch()
    # is it a valid path?
    p.check()

    # graph the path
    p.plot(idx)





if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
