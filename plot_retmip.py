#!/usr/bin/env python3
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


# data import
data = np.genfromtxt('retmip_false.out')
depth = data[0,:][data[0,:] > -9999.0]
dens  = data[1,:][data[1,:] > -9999.0]
temp  = data[2,:][data[2,:] > -9999.0]


data = np.genfromtxt('../RetMIP/RetMIP_Dye-2_long_init.dat')
inp_depth = data[:,0]
inp_dens  = data[:,1]
inp_temp  = data[:,2]


# plot
fig = plt.figure()

# density
ax = fig.add_subplot(121)
ax.plot(
    dens, depth,
    color='#5e81ac', alpha=0.7,
)
ax.plot(
    inp_dens, inp_depth,
    color='#bf616a', alpha=0.7,
)

# temperature
ax = fig.add_subplot(122)
ax.plot(
    temp, depth,
    color='#5e81ac', alpha=0.7,
)
ax.plot(
    inp_temp, inp_depth,
    color='#bf616a', alpha=0.7,
)

# things about the plot
plt.tight_layout()
plt.show()
plt.close()
