#! /usr/bin/env python

"""
File: generate_scatterplot.py
Copyright (c) 2016 Michael Seaman
License: MIT

Description: For use in pt. 3 - takes the maximum of 
a range of x values and plots them with their value of c
"""

import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from final import findmaxima

fig, ax = plt.subplots(nrows = 1, ncols = 1)

c_values = np.arange(2, 6, .01)
for c in c_values:
    sys.stdout.write('.')
    sys.stdout.flush()
    maxima = findmaxima(c)
    for m in maxima:
        ax.plot(c, m, '.')
ax.set_xlabel("c")
ax.set_ylabel("x maxima")
ax.set_ylim((3, 12))
ax.set_title("Maxima on x vs c")
ax.grid(True)
fig.savefig("Scatterplot.png")
plt.close(fig)