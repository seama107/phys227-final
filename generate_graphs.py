#! /usr/bin/env python

"""
File: generate_graphs.py
Copyright (c) 2016 Michael Seaman
License: MIT

Description: For final exam - simply generates the needed graphs
"""

from final import Rossler

import sys
import subprocess
import numpy as np

# c_values = [2, 3, 4, 4.15, 4.2, 5.7]
# for c in c_values:
#     Rc_2 = Rossler(c)
#     print "Running RK4 for c = %f ..." % c
#     Rc_2.run()
#     print "Plotting ..."
#     Rc_2.plotx()
#     Rc_2.plotxy()
#     Rc_2.plotxyz()



# For generating a gif

c_values = np.arange(2, 8, .1)
for c in c_values:
    Rc_2 = Rossler(c)
    print "Running RK4 for c = %f ..." % c
    Rc_2.run()
    print "Plotting ..."
    Rc_2.plotxyz()

subprocess.call("convert -delay 10 -loop 0 *xyz.png Rossler.gif", shell=True)
#subprocess.call("rm *.png", shell=True)

print "Complete!"
