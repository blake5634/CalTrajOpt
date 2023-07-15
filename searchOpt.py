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

#import matplotlib.pyplot as plt

from pympler.classtracker import ClassTracker

##class MemStats(Stats):
#class MemStats:
    #def __init__(self):
        ##Stats.__init__(self)
        #classes = [cto.Cm,cto.point3D,cto.path3D,bd.datafile]
        #self.tracker = ClassTracker()
        #for c in classes:
            #self.tracker.track_class(c)
    #def do_append(self, label):
        #self.tracker.create_snapshot(label)
    #def do_print(self):
        #stats = self.tracker.stats
        #stats.print_summary()
        #stats.dump_stats('pympler.stats')

## example 2
#tracker = ClassTracker()
#tracker.track_object(factory)
#tracker.track_class(Employee)
#tracker.create_snapshot()
#populate_factory(factory)
#tracker.create_snapshot()
#tracker.stats.print_summary()


tracker = ClassTracker()
def start_tracker(classes):
    for c in classes:
        tracker.track_class(c)

def mem_snap(str):
    print(' ... click ...')
    tracker.create_snapshot(str)
    print('Memsnap: ',str)
def mem_report():
    tracker.stats.print_summary()

def main(args):

    #prof
    start_tracker([cto.Cm, cto.point3D,cto.path3D,bd.datafile])

    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    c1 = cto.Cm()  # cost matrix

    #prof
    mem_snap('create Cm')

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


    #prof
    mem_snap('Cm is filled')

    #
    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D(c1)



    q = input('Enter your research question: ')

    #prof
    mem_snap('created path3D()')

    # create the datafile:
    df = bd.datafile('6Dsearching','BH','simulation')
    df.set_folders('/home/blake/Sync/Research/CalTrajOpt_RESULTS','')

    df.metadata.d['Research Question'] = q

    #prof
    mem_snap('created datafile')
    ##########################################################
    #

    #searchType = 'multi heuristic'
    searchType = 'sampling search'

    nsamp = cto.N**6  # at least one for each starting point(!)
    nsamp = 10 # 1M

    #
    ###########################################################


    p.search(searchType,dfile=df,nsamples=nsamp)

    #prof
    mem_snap('complete the search')

    # is it a valid path?
    #p.check()

    #prof
    mem_report()


    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    #p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
