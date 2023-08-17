
import socket
import json
import math
import numpy as np
import matplotlib.pyplot as plt
import brl_data.brl_data as bd
import datetime
import sys
import os
import random


from pympler.classtracker import ClassTracker

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

class beast:
    def __init__(self):
        self.lst = []

    def add(self,x):
        self.lst.append(x)

N = 1000
M = 5

def main(args):
    tstbst = beast()
    start_tracker([beast])
    for j in range(M):
        tstbst.add([5,6,7,8,9]*12)
        for i in range(len(tstbst.lst)):
            print(tstbst.lst[i])
        mem_snap('added N=1000, ones')
    mem_report()

if __name__ ==  '__main__':
    print('main starting:')

    main(sys.argv)
