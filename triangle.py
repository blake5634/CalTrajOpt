#!/usr/bin/python3

#
#    find violations of the triangle inequality in phase space
#            (IF they exist)
#

import random
import c2to as cto
import sys
import brl_data.brl_data as bd
#import matplotlib.pyplot as plt


#def main2D(args):
    #epsilon = 0.02
    #cto.configure()
    #N=4
    #for it in range(100):
        #i = random.randint(0,N-1)
        #j = random.randint(0,N-1)
        #p1 = cto.point2D(i,j)
        #i = random.randint(0,N-1)
        #j = random.randint(0,N-1)
        #p2 = cto.point2D(i,j)
        #i = random.randint(0,N-1)
        #j = random.randint(0,N-1)
        #p3 = cto.point2D(i,j)

        #t12 = cto.trajectory2D(p1,p2)
        #t23 = cto.trajectory2D(p2,p3)
        #t13 = cto.trajectory2D(p1,p3)

        #ts = [t12,t23,t13]
        #for t in ts:
            #dt = 1.0 # initial
            #t.compute(dt) #get coeffs
            #t.constrain_A()    #solve dt to constrain for Amax
            #a = t.timeEvolution(ACC_ONLY=True)
            ##compute traj costs
            #t.cost_e(a)
            #t.cost_t()
        #if not t12.computed:
            #print('somethings wrong')
            #quit()

        #v1 = abs(t12.t_cost+t23.t_cost-t13.t_cost)/t12.t_cost
        #if v1 > epsilon and (t12.t_cost+t23.t_cost) < t12.tcost :
            #print(it,': time cost violated triangle inequality:')
            #print('    p1:',p1)
            #print('    p2:',p2)
            #print('    p3:',p3)
        #v2 = abs( t12.e_cost+t23.e_cost - t13.e_cost)/t12.e_cost
        #if v2>epsilon and (t12.e_cost+t23.e_cost) < t13.e_cost:
            #print(it,': energy cost violated triangle inequality:')
            #print('    p1:',p1)
            #print('    p2:',p2)
            #print('    p3:',p3)

def main3D(args):
    epsilon = 0.0002
    cto.configure()
    N=4
    for it in range(100000):
        i = random.randint(0,N**6-1)
        p1 = cto.point3D(cto.getcoord(i))
        i = random.randint(0,N**6-1)
        j = random.randint(0,N-1)
        p2 = cto.point3D(cto.getcoord(i))
        i = random.randint(0,N**6-1)
        p3 = cto.point3D(cto.getcoord(i))

        t12 = cto.trajectory3D(p1,p2)
        t23 = cto.trajectory3D(p2,p3)
        t13 = cto.trajectory3D(p1,p3)

        ts = [t12,t23,t13]
        for t in ts:
            t.constrain_A()    #solve dt to constrain for Amax
            a = t.timeEvolution(ACC_ONLY=True)
            #compute traj costs
            t.cost_e(a)
            t.cost_t(a)
        if not t12.computed:
            print('somethings wrong')
            quit()

        v1 = abs(t12.t_cost+t23.t_cost-t13.t_cost)/t13.t_cost
        if v1 > epsilon and (t12.t_cost+t23.t_cost) < t13.t_cost :
            print(it,': time cost violated triangle inequality:')
            print('    p1:',p1)
            print('    p2:',p2)
            print('    p3:',p3)
        #v2 = abs( t12.e_cost+t23.e_cost - t13.e_cost)/t13.e_cost
        #if v2 > epsilon and (t12.e_cost+t23.e_cost) < t13.e_cost:
            #print(it,': energy cost violated triangle inequality:')
            #print('    p1:',p1)
            #print('    p2:',p2)
            #print('    p3:',p3)

if __name__ ==  '__main__':
    print('main starting:')
    main3D(sys.argv)
    #main2D(sys.argv)
