#!/usr/bin/python3

import logging
import c2to as cto
# these are needed for unpickling
#from c2to import Cm as Cm
#from c2to import trajectory3D as trajectory3D
#from c2to import point3D as point3D
import sys
import os
import pickle
import brl_data.brl_data as bd
import datetime as dt

#import matplotlib.pyplot as plt

from pympler.classtracker import ClassTracker
from pympler import asizeof

tracker = ClassTracker()
def start_tracker(classes):
    pass
    # commenting these out for speed
    #for c in classes:
        #tracker.track_class(c,resolution_level=2)

def mem_snap(str):
    print(' ... click ...')
    #tracker.create_snapshot(str)  # comment out for speed
    #print('Memsnap: ',str)

def mem_report():
    pass
    #tracker.stats.print_summary()

def main(args):


    q = input('Enter your research question: ')

    # create the datafile:
    df = bd.datafile('6Dsearching','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')

    df.metadata.d['Research Question'] = q


    # memory profiling
    mem_snap('created datafile')


    # memory profiling
    #start_tracker([cto.Cm, cto.point3D,cto.path3D,bd.datafile])
    start_tracker([cto.search_from_curr_pt,cto.path3D])

    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    pts = cto.setupPoints()   # just store points instead of cost matrix Cm
    #c1 = cto.Cm(df = df)  # cost matrix




    # memory profiling
    mem_snap('created path3D()')


    #
    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D()

    ##########################################################
    #
    #  Grid Type:
    RANDOMGRID = True   # remove the pickle file when changed!!!

    #
    # 1D only
    #searchtype = 'brute force'
    # 4x4 and 6D:
    #searchType = 'multi heuristic'
    searchType = 'sampling search'

    # approx 1M samples
    nsamp = 4*cto.N**6  # at least one for each starting point(!)
    #nsamp = N**6 # 1M

    #nsamp = 60 # for testing

    # cost:
    cts = ['time','energy']
    cto.costtype = 'energy'
    #cto.costtype = 'time'
    #
    ###########################################################
    if RANDOMGRID:
        df.metadata.d['Random Grid'] = True
    else:
        df.metadata.d['Random Grid'] = False

    # memory profiling
    ## memory profiling
    #mem_snap('populate Cm')

    #pname = 'c1Costs.pickle'
    #if os.path.exists(pname):
        #print('loading precomputed cost matrix   ...')
        #f = open(pname, 'rb')
        #c1 = pickle.load(f)
        #if RANDOMGRID:
            #c1.set_GridRandomize()  # needed? redundant??
        #print(' loading completed')
    #else:
        #print('no stored data: computing cost matrix')
        #if RANDOMGRID:
            ##print('    ... bug: always have to recompute in Random Grid mode!!!')
            #c1.set_GridRandomize()  # random pts instead of grid
        #c1.fill()
        #f = open(pname,'wb')
        #pickle.dump(c1,f)


    #mem_snap('Cm is filled')



    p.search(searchType,dfile=df,nsamples=nsamp,profiler=mem_snap)
    # search will close the datafile

    # memory profiling
    mem_snap('completed the search')

    # is it a valid path?
    #p.check()

    # memory profiling
    mem_report()

    ####  keep a "log book"
    dim = '6D'
    logdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    notes = '{:}, n={:}, {:}, cost: {:4.2f} ({:})'.format(searchType, nsamp,dim, df.metadata.d['Min Cost'], cto.costtype)
    now = dt.datetime.now()
    dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
    logentry = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, df.metadata.d['Research Question'])
    f =open(logdir+'work_logbook.txt','a')
    print(logentry,file=f)
    f.close()
    ####

    print('Completed: see results at ',df.hashcode)

    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    #p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
