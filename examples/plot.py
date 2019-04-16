import matplotlib.pyplot as plt
from math import *
import os
import sys

files = []
fpoints = []
if len(sys.argv) > 1:
    files = sys.argv[1:]
else:
    files = ['-']

fig, ax = plt.subplots(nrows=len(files), ncols=1, gridspec_kw={'hspace':1}, figsize=(8, 6), dpi=160)
if len(files) == 1:
    ax = [ ax ]

for idx, fname in enumerate(files):
    title = fname 
    opts = ''
    if fname.count(':') == 1:
        fname, title = fname.split(':')
    elif fname.count(':') == 2:
        fname, title, opts = fname.split(':')


    fl = sys.stdin if fname == '-' else open(fname, 'r')
    points = []
    for line in fl:
        points += [ (*map(float, line.split(' ')), ) ]
    fpoints.append(points)

    if opts.find('R') != -1:
        opts,ropts = opts.split('R')
        ax[idx].plot(*zip(*fpoints[idx-1]), ropts)

    ax[idx].set_title(title)
    X, Y = zip(*points)
    ax[idx].plot(X,Y, opts)

    fl.close()

if 'GUI' in os.environ:
    plt.show(block=True)

fig.savefig('fig.png') 
