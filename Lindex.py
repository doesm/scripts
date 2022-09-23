#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 11:50:29 2019

@author: Marcelo Caparotta

Lindemann Index.
"""

import argparse
import subprocess
import pandas as pd
import math

parser = argparse.ArgumentParser()
parser.add_argument('trajectory', help='XTC Trajectory')
parser.add_argument('start', type=int, help='First atom')
parser.add_argument('end', type=int, help='Last atom')
args = parser.parse_args()

dv = ""
dp = ""
with open('lindex.dat', 'w') as f:
    for i in range(args.start, args.end+1):
        for j in range(args.start, args.end+1):
            if j != i:
                f.write('d' + str(i)+'-'+str(j)+': DISTANCE ATOMS=' +
                        str(i) + ',' + str(j) + '\n')
                f.write('p' + str(i)+'-'+str(j)+': CUSTOM ARG=d' +
                        str(i)+'-'+str(j)+' FUNC=x*x PERIODIC=NO\n')
                if j == args.start:
                    dv = 'd'+str(i)+'-'+str(j)
                    dp = 'p'+str(i)+'-'+str(j)
                else:
                    dv = dv+',d'+str(i)+'-'+str(j)
                    dp = dp+',p'+str(i)+'-'+str(j)
        f.write('PRINT ARG='+dv+' FILE=dv'+str(i)+'\n')
        f.write('PRINT ARG='+dp+' FILE=dp'+str(i)+'\n')
f.closed

subprocess.run(["plumed", "driver", "--plumed",
               "lindex.dat", "--mf_xtc", args.trajectory])

sumj = 0
sumqi = 0
N = args.end+1-args.start
for i in range(args.start, args.end + 1):
    # print(i)
    df = pd.read_csv('dv'+str(i)+'', header=None, skiprows=1,
                     skipinitialspace=True, index_col=0, sep='\s+', dtype=float)
    pf = pd.read_csv('dp'+str(i)+'', header=None, skiprows=1,
                     skipinitialspace=True, index_col=0, sep='\s+', dtype=float)
    for j in range(1, N):
        # print(j)
        sumj = sumj+math.sqrt(pf.mean(axis=0)
                              [j]-df.mean(axis=0)[j]**2)/df.mean(axis=0)[j]
    qi = (1/(N-1))*sumj
    print("Atom ", i, " = ", qi)
    sumj = 0
    sumqi = sumqi+qi

q = 1/N*sumqi
print("System-averaged Lindemann index: ", q)
