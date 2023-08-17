#!/usr/bin/python3

import glob
import os
import c2to as cto
import sys
import brl_data.brl_data as bd
import re
import datetime as dt

hashmatcher = r'[^a-f,^0-9]*([a-f, 0-9]{8})[^a-f,^0-9]*' # 8character hex hash
hmc = re.compile(hashmatcher)

def main(args):

    #  all the dirs we might find data files
    dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D-Round2',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom',
            '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing'
            ]
    flags = ['nothing']
    if len(args) > 1: # cmd line delete specific hash files and purge from log
        hlist = args[1:]
        #print('Test0: hlist: ',hlist)
        if hlist[0] == '-listOnly':
            #print(' I caught listOnly flag!')
            flags += ['listOnly']
            hlist = hlist[1:]  # cut first entry
            #print('Test1: hlist: ',hlist)
        hrem = purgeFilesbyHash(hlist,dirs,flags=flags)
        for h in hrem:
            log_record_deletions(h,'file')
        hrem = purgeLogsbyHash(hlist,flags=flags)
        for h in hrem:
            log_record_deletions(h,'log entry')
    else:
        # without cmd line args:
        #   auto mode based on metadata and Research Question
        datadir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS'
        hashs = cleanOutDatafiles([datadir])
        if len(hashs) > 0:
            print('          Now cleaning removed files from logs:')
            hrem = purgeLogsbyHash(hashs)
            for h in hrem:
                log_record_deletions(h,'log entry')
        else:
            print('          No log file cleaning is needed. Logfile unmodified')
    quit()

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
        FILESWEREREMOVED = False
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
                    hashesRemoved.add(getHashFromFilename(fn))
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
                FILESWEREREMOVED = True
                for fn in remlistCSV:
                    print(' ... removing ',fn)
                    hashesRemoved.add(getHashFromFilename(fn))
                    os.remove(fn)

            else:
                print('removing canceled')
        if FILESWEREREMOVED:
            print('\n          Done removing files\n')
        else:
            print('\n          No files were matched or removed\n')
        return list(hashesRemoved)

def purgeFilesbyHash(hlist,dirs,flags=['None']):
    hashesRemoved = set()
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
    if 'listOnly' in flags:
        return []
    x = input('\n\n          OK to delete these files?? ... (y/N)')
    if x.lower() == 'y':
        for fd in fds:
            fn = fd[1]
            print(' ... removing ',fn)
            hashesRemoved.add(getHashFromFilename(fn))
            os.remove(fd[0]+'/'+fd[1])
        return list(hashesRemoved)
    else:
        print('removing canceled')
        return []


def getHashFromFilename(fn):
        result = hmc.findall(fn)
        if len(result) < 1:
            print('getHashFromFilename result:',result)
            print(f'No hash found in filename: {fn} (should have {hver1})')
            quit()
        if len(result) > 1:
            print(f'getHashfromFilename - warning: {len(result)} hashes found only first one used')
        hashResult = result[0]
        #print(f'result: {result}, hv2: {hashResult}')
        return hashResult


def purgeLogsbyHash(hlist,flags=['None']):
    print(f'\nChecking  {hlist} in log files.')
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
        if 'listOnly' in flags:
            pass
            #print(f'\n    Matched {len(deletelines)} {hlist} entries from {justname}')
        else:
            print(f'\n    Deleting {len(deletelines)} {hlist} entries from {justname}')
        f.close()

        # explain what you will do
        print(f'{justname}:')
        i=0
        for l in deletelines:
            i+=1
            l = l.strip()
            print(f'    {i}: [{l}]')

        # confirm y/N ??
        if 'listOnly' not in flags:
            x = input(f'\n   OK to delete {len(deletelines)} entries from {justname}? (y/N): ')
            if x.lower() == 'y':
                # now modify the file
                f = open(lf,'w')
                for l in keeplines:
                    print(l.strip(),file=f)
                f.close()
                print(f' {len(deletelines)} lines purged from {justname}.')
            else:
                print(' purging logs canceled.')
    if 'listOnly' in flags:
        print('\n\n             List only, no actions were taken.')
    return [] # we don't log hashes deleted from log files

def log_record_deletions(hash, dtype):
    logfile = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/deletions.txt'
    f = open(logfile, 'a')
    now = dt.datetime.now()
    dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
    logline = f'{dtstring}: deleted {dtype} for {hash}'
    print(logline,file=f)
    f.close()

def logentry(df,notes):
    logdir = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/'
    logfilename = 'work_logbook.txt'
    q = df.metadata.d['Research Question']
    if len(q)>0 and 'debug' not in q:  # skip junk files
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        logline = '{:}, {:}, {:},  RQ: {:}'.format(dtstring,df.hashcode, notes, q)
        with open(logdir+logfilename,'a') as f:
            print(logline,file=f)
            f.close()
            print('added log entry to: ',logdir+'work_logbook.txt')
    else:
        print(f'debugging detected. {df.hashcode} will not be logged to {logdir+logfilename}')



if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
