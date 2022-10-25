#!/usr/bin/env python3
#-*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt


comma2dot = lambda x: float(x.replace(',', '.'))



# kohnen data
data = './kohnen04init.dat'
data = np.genfromtxt(data)
kohnen_depth = data[:,0]
kohnen_dens  = data[:,1]


# simulation results
data = './kohnen04.out'
data = np.genfromtxt(data, skip_footer=1)
depth = data[:,0]
dens  = data[:,1]
temp  = data[:,2]
depth = depth - depth[-1]

# bagmean
dz = 1.0
mean_depth = np.array([])
mean_dens  = np.array([])
for z in np.arange(depth[-1], depth[0], -dz):
    zmin = z - (0.5 * dz)
    zmax = z + (0.5 * dz)

    mask = np.logical_and(
        (depth >= zmin),
        (depth <  zmax),
    )

    mean_depth = np.append(mean_depth, z)
    mean_dens  = np.append(mean_dens, np.mean(dens[mask]))



# test plot
fig = plt.figure(figsize=(7,9))
ax = fig.add_subplot(111)

ax.plot(
    kohnen_dens, kohnen_depth,
    color='#5e81ac',
    alpha=0.7, linewidth=3.0,
    label='Kohnen (init)'
)
ax.plot(
    dens, depth,
    #dens[depth >= -25.0], depth[depth >= -25.0],
    color='#bf616a',
    alpha=0.7, linewidth=3.0,
    label='simulation result',
)
#ax.plot(
#    mean_dens, mean_depth,
#    #mean_dens[mean_depth >= -25.0], mean_depth[mean_depth >= -25.0],
#    color='#bf616a',
#    alpha=0.7, linewidth=2.0
#)

ax.spines['bottom'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel(r'density $(\mathrm{kg \, m}^{-3})$')
ax.set_ylabel(r'depth $(\mathrm{m})$')
ax.xaxis.tick_top()
ax.xaxis.set_label_position('top')
ax.legend(loc=1, frameon=False)

# things about the plot
plt.tight_layout()
plt.savefig('kohnen04long_density.png', dpi=300)
plt.show()
plt.close()
