#!/usr/bin/python3

import logging
import c2to as cto
# these are needed for unpickling
from c2to import Cm as Cm
from c2to import trajectory3D as trajectory3D
from c2to import point3D as point3D
import sys
import os
import pickle
import brl_data.brl_data as bd

#import matplotlib.pyplot as plt


def main(args):
    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    c1 = Cm()  # cost matrix


    pname = 'c1Costs.pickle'
    if os.path.exists(pname):
        print('loading precomputed cost matrix   ...')
        f = open(pname, 'rb')
        c1 = pickle.load(f)
        print(' loading completed')
    else:
        print('no stored data: computing cost matrix')
        c1.fill()
        f = open(pname,'wb')
        pickle.dump(c1,f)


    #
    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D(c1)

    # create the datafile:
    df = bd.datafile('6Dsearching','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')

    q = input('Enter your research question: ')
    df.metadata.d['Research Question'] = q

    ##########################################################
    #

    #searchType = 'multi heuristic'
    searchType = 'sampling search'

    nsamp = cto.N**6  # at least one for each starting point(!)
    nsamp = 1000000  # 1M

    #
    ###########################################################


    p.search(searchType,dfile=df,nsamples=nsamp)

    # is it a valid path?
    #p.check()

    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    #p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
