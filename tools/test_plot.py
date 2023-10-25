#!/usr/bin/env python3
#-*- coding: utf-8 -*-


import glob
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


DATA_DIR = '../tfm_output/*'
FORCING_FILE = '../RetMIPexamples/FA_forcing.dat'


def importData(data_file):
    
    data = np.genfromtxt(data_file)

    depth        = data[:,0]
    density      = data[:,1]
    temperature  = data[:,2]
    grain_radius = data[:,3]
    age          = data[:,4]
    liquid_water = data[:,5]
    
    #liquid_water = liquid_water / ((1.0 - (density / 917.0)) * (917.0 / 1000.0))

    return {
        'depth'        : depth,
        'density'      : density,
        'temperature'  : temperature,
        'grain_radius' : grain_radius,
        'age'          : age,
        'liquid_water' : liquid_water,
    }


def importForcing(forcing_file):
    
    data = np.genfromtxt(forcing_file, delimiter=',')
    model_time        = data[:,1] / (3600.0 * 24.0)
    surf_temperature  = data[:,2]
    surf_density      = data[:,3]
    surf_solid_accu   = data[:,4] * (3600.0 * 24.0 * 365.0)
    surf_liquid_accu  = data[:,5] * (3600.0 * 24.0 * 365.0)
    surf_grain_radius = data[:,6]

    return {
        'model_time'        : model_time,
        'surf_temperature'  : surf_temperature,
        'surf_density'      : surf_density,
        'surf_solid_accu'   : surf_solid_accu,
        'surf_liquid_accu'  : surf_liquid_accu,
        'surf_grain_radius' : surf_grain_radius,
    }


if __name__ == '__main__':
    
    # forcing data import
    forcing = importForcing(FORCING_FILE)
    min_max_temp = (
        min(forcing['surf_temperature']),
        max(forcing['surf_temperature']),
    )
    min_max_solid_accu = (
        min(forcing['surf_solid_accu']),
        max(forcing['surf_solid_accu']),
    )
    min_max_liquid_accu = (
        min(forcing['surf_liquid_accu']),
        max(forcing['surf_liquid_accu']),
    )
    min_max_accu = (
        min((min_max_solid_accu[0], min_max_liquid_accu[0])),
        max((min_max_solid_accu[1], min_max_liquid_accu[1])),
    )


    t = 0
    # loop over all output files
    for data_file in sorted(glob.glob(DATA_DIR)):
        
        # data import
        data = importData(data_file)

        mask = np.logical_and(
            (forcing['model_time'] >= (forcing['model_time'][t] - 10)),
            (forcing['model_time'] <= (forcing['model_time'][t] + 10)),
        )
        mask = np.ones_like(forcing['model_time'], dtype=bool)


        fig = plt.figure(figsize=(16,9))
        gs = gridspec.GridSpec(3, 3, height_ratios=[1,1,5])

        # density
        ax = fig.add_subplot(gs[2,0])
        ax.plot(
            data['density'], data['depth'],
            color='#3b4252', linewidth=2.0,
        )
        ax.set_xlabel(r'density $(\mathrm{kg \, m}^{-3})$')
        ax.set_ylabel(r'depth $(\mathrm{m})$')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_xlim(250.0, 920.0)

        # liquid water
        ax = fig.add_subplot(gs[2,1])
        ax.plot(
            data['liquid_water'], data['depth'],
            color='#3b4252', linewidth=2.0,
        )
        #ax.set_xlabel(r'effective liquid water content $(\mathrm{1})$')
        ax.set_xlabel(r'liquid water content $(\mathrm{1})$')
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([])
        ax.set_yticklabels([])
        #ax.set_xlim(-0.05, 1.05)

        # temperature
        ax = fig.add_subplot(gs[2,2])
        ax.plot(
            data['temperature'], data['depth'],
            color='#3b4252', linewidth=2.0,
        )
        ax.axvline(273.15, color='#3b4252', linewidth=1.5, linestyle='--')
        ax.set_xlabel(r'temperature $(\mathrm{K})$')
        ax.set_ylabel(r'depth $(\mathrm{m})$')
        ax.yaxis.tick_right()
        ax.yaxis.set_label_position('right')
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        #ax.set_xlim(180.0, 275.0)



        # temperature
        ax = fig.add_subplot(gs[0,:])
        ax.plot(
            forcing['model_time'][mask], forcing['surf_temperature'][mask],
            color='#8fbcbb', linewidth=2.0,
        )
        ax.axvline(forcing['model_time'][t], color='#bf616a', linewidth=2.0)

        ax.set_ylabel(r'surface' + '\n' + 'temperature $(\mathrm{K})$')
        ax.set_xticks([])
        ax.set_xticklabels([])
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(min_max_temp[0], min_max_temp[1])

        # accumuation / ablation
        ax = fig.add_subplot(gs[1,:])
        ax.plot(
            forcing['model_time'][mask], forcing['surf_solid_accu'][mask],
            color='#5e81ad', linewidth=2.0,
        )
        ax.plot(
            forcing['model_time'][mask], forcing['surf_liquid_accu'][mask],
            color='#bf616a', linewidth=2.0,
        )
        ax.axvline(forcing['model_time'][t], color='#bf616a', linewidth=2.0)

        ax.set_xlabel(r'model time $(\mathrm{d})$')
        ax.set_ylabel(r'accumulation' + '\n' + r'$(\mathrm{m\,weq.\,a}^{-1})$')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_ylim(min_max_accu[0], min_max_accu[1])

        # things about the plot
        plt.savefig('./plots/plot%06i.png' % t, dpi=100)
        #plt.show()
        plt.close()

        t += 1
