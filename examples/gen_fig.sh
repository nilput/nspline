#!/bin/bash
make && \
./test_prog 200 0 10  < vals.txt > interp.txt && \
/usr/bin/env python3 plot.py 'vals.txt:Datapoints:-o' 'interp.txt:Cubic spline interpolation:bRob'
