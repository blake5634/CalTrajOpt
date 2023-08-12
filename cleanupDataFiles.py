#!/usr/bin/python3

import glob
import os
import c2to as cto
import sys
import brl_data.brl_data as bd

def main(args):

    #  all the dirs we might find data files
    dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D-Round2',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing'
            ]
    if len(args) > 1: # cmd line delete specific hash files and purge from log
        hlist = args[1:]
        purgeFilesbyHash(hlist,dirs)
        purgeLogsbyHash(hlist)
        quit()

    # without cmd line args:
    #   auto mode based on metadata and Research Question
    datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
    hashs = cleanoutDatafiles([datadir])
    print(' Now cleaning removed files from logs:')
    purgeLogsbyHash(hashs)

def cleanOutDatafiles(dirs):  # based on metadata
    #####################################################################
    #
    # get all the files
    #
    # thanks stackoverflow
    for datadir in dirs:
        hashesRemoved = set()
        files = list(glob.glob(datadir + '/' + "*"))
        files.sort(key=lambda x: os.path.getmtime(x),reverse=True) # newest first
        filenameroots = []
        if len(files) ==0:
            cto.error('No brl_data files found')
        for f in files:
            #print('found: ',f)
            if '.csv' in f :
                filenameroots.append(str(f).replace('.csv',''))
            if  '_meta.json' in f:
                filenameroots.append(str(f).replace('_meta.json',''))
        # dict.fromkeys better than set because preserves order (thanks google)
        filenameroots = list(dict.fromkeys(filenameroots)) # elim dupes


        #
        #  OK - now lets check all metadatas
        #
        remlist = []
        remlistCSV = []
        for fnr in filenameroots:
            mdname = fnr+'_meta.json'
            if(os.path.exists(mdname)):  # if we can find a metadata file
                df = bd.datafile('', '','')  # open it with blank title info
                df.set_folders(datadir,'')        # set these to wherever you want to open datafiles
                df.open(mode='r',tname=mdname)
                #df.metadata.polish()  # convert metadata from strings to useful types
                df.close()
                # RQ key came in two forms: with and without a space
                try:
                    q = df.metadata.d['Research Question']
                    if len(q)==0 or 'debug' in q.lower(): # there is no real RQ
                        remlist.append(fnr+'_meta.json')
                        remlist.append(fnr+'.csv')
                except:
                    pass
                try:
                    q = df.metadata.d['ResearchQuestion']
                    if len(q)==0 or 'debug' in q.lower(): # there is no real RQ
                        remlist.append(fnr+'_meta.json')
                        remlist.append(fnr+'.csv')
                except:
                    pass
            else:
                print('non existent metadata file: ',fnr)
                if 'ties' not in fnr:  # we have special csv files for ties data
                    remlistCSV.append(fnr+'.csv')

        if len(remlist) > 0:
            print('\n\nPlanning to remove the following file names (.jason and .csv):')
            print('     Reason: Research Question is empty or contains "debug"')
            for fn in remlist:
                print(fn)

            x = input('\n\n          OK to delete these files?? ... (y/N)')
            if x.lower() == 'y':
                for fn in remlist:
                    # keep around a set of the hashes removed
                    hashsRemoved.add(fn.split('_')[-4])   # fragile!!
                    print(' ... removing ',fn)
                    os.remove(fn)
            else:
                print('removing canceled')

        if len(remlistCSV) > 0:
            print('\n\nPlanning to remove the following .csv files')
            print('     Reason: they have no metadata file (often ^C or crash)')

            for fn in remlistCSV:
                print(fn)

            x = input('\n\n          OK to delete these files?? ... (y/N)')

            if x.lower() == 'y':
                for fn in remlistCSV:
                    print(' ... removing ',fn)
                    os.remove(fn)
            else:
                print('removing canceled')
        print('Done removing files')
        return list(hashsRemoved)

def purgeFilesbyHash(hlist,dirs):
    myf = bd.finder()
    myf.set_dirs(dirs)
    fds = []
    for h in hlist:
        keys = [h]
        fds += myf.findh(keys) # get all files with this hash
    print(f' I found {len(fds)} files matching {hlist}:')
    i=0
    for fd in fds:
        i+=1
        print(f'    {i}: [{fd[0]}/{fd[1]}]')
    x = input('\n\n          OK to delete these files?? ... (y/N)')
    if x.lower() == 'y':
        for fd in fds:
            print(' ... removing ',fd[1])
            os.remove(fd[0]+'/'+fd[1])
    else:
        print('removing canceled')
    #

def purgeLogsbyHash(hlist):
    x = input(f'\n\n          OK to purge {hlist} from log files?? ... (y/N)')
    if x.lower() == 'y':
        # now purge entries from the logs
        wlog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/work_logbook.txt'
        ilog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/image_log.txt'
        logs = [wlog, ilog]
        for lf in logs:
            justname = lf.split('/')[-1]
            f = open(lf,'r')
            keeplines = []
            deletelines = []
            for line in f:
                save = True
                for h in hlist:
                    if h in line:
                        save = False
                if save:
                    keeplines.append(line)
                else:
                    deletelines.append(line)
            print(f'\n Deleting {len(deletelines)} {hlist} entries from {justname}')
            f.close()
            for l in deletelines:
                print('deleting: ',l.strip())

            # confirm y/N ??

            # now do it
            f = open(lf,'w')
            for l in keeplines:
                print(l.strip(),file=f)
            f.close()
            print(f' {len(deletelines)} lines purged from {justname}.')
    else:
        print(' purging logs canceled.')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
