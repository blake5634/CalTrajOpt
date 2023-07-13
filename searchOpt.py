#!/usr/bin/python3

import logging
import c2to as cto
import sys
import brl_data.brl_data as bd

#import matplotlib.pyplot as plt


def main(args):
    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    Cm = Cm()  # cost matrix
    Cm.fill()  # fill in with trajectories and costs

    #
    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D(Cm, adv=True)

    # create the datafile:
    df = bd.datafile('6Dsearching','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')


    ##########################################################
    #

    searchType = 'multi heuristic'
    #searchType = 'sampling search'

    nsamp = cto.N**6  # at least one for each starting point(!)

    #
    ###########################################################


    p.search(searchType,dfile=df,nsamples=nsamp)

    # is it a valid path?
    #p.check()

    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
