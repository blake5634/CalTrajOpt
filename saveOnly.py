#!/usr/bin/python3

import glob
import os
import c2to as cto
import sys
import brl_data.brl_data as bd
import datetime as dt
import shutil as shu

DEBUG = False
###ACTION = 'delete'   # actually delete old files
ACTION = 'archive'  # save old files to trash
#ACTION = 'simulate'  # just simulate actions, don't do anything


def main(args):

    ##  all the dirs we might find data files
    #dirs = ['/home/blake/Sync/Research/CalTrajOpt_RESULTS',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D-Round2',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/1D_data/Gold',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/PointSetsRandom',
            #'/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing'
            #]

    #actionlogfile = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/deletions.txt'

    #wlog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/work_logbook.txt'
    #ilog = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/writing/image_log.txt'
    #logs = [wlog, ilog]
    #trashcan = '/home/blake/Sync/Research/CalTrajOpt_RESULTS/Archive'

    if True:
        # for testing
        dirs = ['/home/blake/Ptmp/CalTrajOpt/testing/folder1',
                '/home/blake/Ptmp/CalTrajOpt/testing/folder2',
                '/home/blake/Ptmp/CalTrajOpt/testing/folder3'
                ]

        logs = ['/home/blake/Ptmp/CalTrajOpt/testing/log1.txt',
                '/home/blake/Ptmp/CalTrajOpt/testing/log2.txt' ]

        actionlogfile = '/home/blake/Ptmp/CalTrajOpt/testing/deletions.txt'
        trashcan = '/home/blake/Ptmp/CalTrajOpt/testing/Trash'

    locinfo = { 'dirs' : dirs,
                'logs' : logs,
                'actionlogfile' : actionlogfile,
                'trashcan' : trashcan }

    flags = ['nothing']

    if len(args) != 2:
        bd.brl_error('usage: >saveOnly.py   URIfile')


    # 0) collect the URLs to save
    uriFile = args[1]
    try:
        f = open(uriFile,'r')
    except:
        bd.brl_error(f'URI file: {uriFile} not found')
    URIs2Save = []
    for l in f:
        l = l.strip()
        if len(l) < 6 or len(l)>8:
            bd.brl_error(f'URI:{l} is less than minimum 6 characters or > 8')
        for c in l:
            if c not in 'abcdef0123456789':
                bd.brl_error(f'illegal URI charcter: {c} in {l}')
        URIs2Save.append(l)
    n = len(URIs2Save)

    print(f' I found {n} URIs to save.  First 3:')
    for i in range(3):
        print(URIs2Save[i])
    x = input('  ...  OK? ... (y/N)')
    if x.lower() != 'y':
        quit()

    # 1) find all files NOT in the URIs2Save list
    fdList = []
    hset = set()
    for d in dirs:
        files = os.listdir(d) # file names
        for fn in files:
            theseHs = bd.getHashFromFilename(fn)
            for thishash in theseHs:
                #print('processing: ', thishash)
                # do not collect the hash if it is "None" or is a keeper
                if thishash not in URIs2Save and thishash is not None:
                    #print('   old file hash list: appending: ',thishash)
                    fdList.append([d,f])
                    hset.add(thishash)

    targetFileHashList = []
    n = 0
    for hl in list(hset):
        n+=1
        targetFileHashList.append(hl)  # these are 1-member lists: [ filename ]

    print(f'We have identified {len(targetFileHashList)} hashes to {ACTION}.')

    # 2) purge those files

    hrem = processFilesbyHash(targetFileHashList,dirs,locinfo)
    #for h in hrem:
        #log_record_deletions(h,ACTION,'file',locinfo)


    # 3) find all log entries NOT having URIs2Save in them (hsetLogLines)

    hsetLogLines = set() # hashes to delete lines from log files
    for fn in logs:
        f = open(fn,'r')
        for line in f:
            hashs = bd.getHashFromFilename(line) # [] if line has no hash in it
            for thishash in hashs:
                if thishash not in URIs2Save and thishash is not None:
                    hsetLogLines.add(thishash)
        f.close()

    # 4) purge the log entries
    hrem = processLogsbyHash(list(hsetLogLines),logs, locinfo, flags=flags)
    quit()

def operation(ACTION, fdir, fname,locinfo):
    sourcefile = fdir+'/'+fname
    if ACTION == 'archive':
        shu.move(sourcefile, locinfo['trashcan'] + '/' + fname)
    elif ACTION == 'delete':
        os.remove(sourcefile)
    elif ACTION == 'simulate':
            pass
    else:
        print('unknown ACTION: ',ACTION)
        quit()

def processFilesbyHash(targets,dirs,locinfo, flags=['None']):
    hashesRemoved = set()
    myf = bd.finder()
    myf.set_dirs(dirs)
    fds = [] # 'fd' is a pair (dir,fname)
    for h in targets:
        keys = [h]
        fds += myf.findh(keys) # get all files with this hash
    print(f' I found {len(fds)} files matching the hash list:')

    i=0
    for fd in fds:  # show the user all the files you will work on
        i+=1
        print(f'    {i}: [{fd[0]}/{fd[1]}]')
    if 'listOnly' in flags:
        return []
    x = input(f'\n\n          OK to {ACTION} these files?? ... (y/N)')
    if x.lower() == 'y':
        for fd in fds:
            fn = fd[1]
            if ACTION == 'delete':
                print(' ... removing ',fn)
            elif ACTION == 'simulate':
                print(' ... simulating ',fn)
            else:
                print(f" ... archiving {fn} to {locinfo['trashcan']} ...")
            for h in bd.getHashFromFilename(fn):
                hashesRemoved.add(h)
            operation(ACTION, fd[0],fd[1],locinfo)
            log_record_deletions(h,ACTION,'file',locinfo['actionlogfile'])
        return list(hashesRemoved)
    else:
        print('removing canceled')
        return []

def processLogsbyHash(hlist,logs,locinfo, flags=['None']):
    print(f'\nChecking  {hlist} in log files.')
    # now purge entries from the logs
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
            if save:  # always keep the lines if simulating
                keeplines.append(line)
            if not save:
                deletelines.append(line)

        if 'listOnly' in flags:
            pass
            #print(f'\n    Matched {len(deletelines)} {hlist} entries from {justname}')
        else:
            # explain what you will do
            verb = f'{ACTION}ing'.replace('ein','in')
            print(f'\n    {verb} {len(deletelines)}  entries from {justname}:')
        f.close()

        i=0
        for l in deletelines:
            i+=1
            l = l.strip()
            print(f'    {i}: [{l}]')

        # confirm y/N ??
        if 'listOnly' not in flags:
            x = input(f'\n   OK to {ACTION} {len(deletelines)} entries from {justname}? (y/N): ')
            if x.lower() == 'y':
                if ACTION != 'simulate':
                    # now modify the file
                    f = open(lf,'w')
                    for l in keeplines:
                        print(l.strip(),file=f)
                    f.close()
                    for l in deletelines:
                        hs = bd.getHashFromFilename(l)  # maybe can be more than one?
                        for h in hs:
                            log_record_deletions(h,ACTION,'log entry',locinfo['actionlogfile'])
                verb = f'{ACTION}ed'.replace('eed','ed')
                print(f' {len(deletelines)} lines {verb} from {justname}.')
            else:
                print(' purging logs canceled.')
    if 'listOnly' in flags:
        print('\n\n             List only, no actions were taken.')
    return [] # we don't log hashes deleted from log files

def log_record_deletions(hash, action, dtype, logfile):
    if action == 'simulate':
        print(f'   simulating recording a {dtype} deletion in {logfile}')
        return
    else:
        f = open(logfile, 'a')
        now = dt.datetime.now()
        dtstring = now.strftime('%Y-%m-%d %H:%M:%S')
        verb = f'{action}ed'.replace('ee','e')
        logline = f'{dtstring}: {verb} {dtype} for {hash}'
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
