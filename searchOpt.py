#!/usr/bin/python3

import logging
import c2to as cto
import sys
#import matplotlib.pyplot as plt


def main(args):
    # create a path and plot it graphically
    cto.configure()
    idx = int(args[1])

    # compute a path:
    #p = cto.path3d(gt,c1)
    p = cto.path3D(adv=True)
    # is it a valid path?
    #p.check()

    # graph the path
    #p.plot(idx)
    # save the path for graphing by animate3D.py
    p.save('3Dsearch')


if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
