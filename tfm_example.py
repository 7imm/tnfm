#!/usr/bin/env python3
#-*- coding: utf-8 -*-


import urllib.request
import numpy as np
import matplotlib.pyplot as plt


# conversation constants
CELSIUS2KELVIN = 273.15
SECONDS_YEAR = (3600.0 * 24.0 * 365.0)

# constant initial data
INIT_GRAIN_RADIUS = 0.0005
INIT_LIQUID_WATER = 0.0

# constant forcing cata
FORCING_DENSITY      = 360.0
FORCING_LIQUID_ACC   = 0.0
FORCING_GRAIN_RADIUS = 0.0005

# model defintions
SOLVE_DENSITY              = 'herron1980'
SOLVE_TEMPERATURE          = 'true'
SOLVE_HEAT_CAPACITY        = 'paterson1994'
SOLVE_THERMAL_CONDUCTIVITY = 'sturm2007'
SOLVE_LIQUID               = 'bucket'
SOLVE_GRAIN_GROWTH         = 'arthern2010'

FORCING_INPUT_FILE         = './tfm_example_forcing.dat'
INITIAL_INPUT_FILE         = './tfm_example_init.dat'

SPINUP                     = -1.0

# Wilhelms, F. (2007). Density of firn core DML94C07_38. Alfred Wegener
# Institue, Helmholtz Centre for Polar and Marine Research, Bremerhaven,
# PANGAEA, https://doi.org/10.1594/PANGAEA.615180
INIT_DENS_URL = 'https://doi.pangaea.de/10.1594/PANGAEA.615180?format=textfile'
INIT_DENS_FILE = 'init_dens.dat'
INIT_DENS_COLS = [0, 1]

# Fernandoy, F., Meyer, H., Oerter, H., Wilhelms, F., Graf, W., and
# Schwander, J. (2010). Borehole temperature profile of ice core
# DML94C07_38. PANGAEA, https://doi.org/10.1594/PANGAEA.753165
INIT_TEMP_URL = 'https://doi.pangaea.de/10.1594/PANGAEA.753165?format=textfile'
INIT_TEMP_FILE = 'init_temp.dat'
INIT_TEMP_COLS = [0, 2]

# Fernandoay, F., Meyer, H., Oerter, H., Wilhelms, F., Graf, W., and
# Schwander, J. (2008). Annual means of d18O and accumulation rates of
# firn core DML09C07_38. Alfred Wegener Institute, Helmholtz Centre for
# Polar and Marine Research, Bremerhaven, PANGAEA,
# https://doi.org/10.1594/PANGAEA.704480
INIT_AGE_URL = 'https://doi.pangaea.de/10.1594/PANGAEA.704480?format=textfile'
INIT_AGE_FILE = 'init_age.dat'
INIT_AGE_COLS = [0, 2]

# Fernandoy, F., Meyer, H., Oerter, H., Wilhelms, F., Graf, W., and
# Schwander, J. (2010). Monthly mean d18O values of ice core DML94C07_38
# from 1981 to 2007. PANGAEA, https://doi.org/10.1594/PANGAEA.753171
FORCING_TEMP_URL = 'https://doi.pangaea.de/10.1594/PANGAEA.753171?format=textfile'
FORCING_TEMP_FILE = 'forcing_temp.dat'
FORCING_TEMP_COLS = [0, 1]

# Fernandoy, F., Meyer, H., Oerter, H., Wilhelms, F., Graf, W., and
# Schwander, J. (2008). Annual means of d18O and accumulation rates of
# firn core DML94C07_38. Alfred Wegener Institute, Helmholtz Centre for
# Polar and Marine Research, Bremerhaven, PANGAEA,
# https://doi.org/10.1594/PANGAEA.704480
FORCING_ACC_URL = 'https://doi.pangaea.de/10.1594/PANGAEA.704480?format=textfile'
FORCING_ACC_FILE = 'forcing_acc.dat'
FORCING_ACC_COLS = [2, 5]



def data_import(file, cols, dtype=float):

    # header identification
    f = open(file, 'r')
    header_length = 0
    for line in f:
        if ('*/' in line):
            header_length += 2
            break

        header_length += 1

    f.close()

    data = np.genfromtxt(
        file,
        skip_header=header_length,
        dtype=dtype,
        usecols=cols,
        delimiter='\t',
    )

    return data



def fix_temperature_dates(dates):
    
    dt = (dates[1] - dates[0])

    for n in range(len(dates) - 1):
        if (round((dates[n+1] - dates[n]), 4) != round(dt, 4)):
            dates[n+1] = dates[n] + dt

    return dates


def generate_config(solve_density, solve_temperature, solve_heat_capacity,
    solve_thermal_conductivity, solve_liquid, solve_grain_growth,
    forcing_input_file, initial_input_file, time_step, spinup):
    
    conf = (
        '&CONFIG\n'
        '  solve_density="%s",\n'
        '  solve_temperature="%s",\n'
        '  solve_heat_capacity="%s",\n'
        '  solve_thermal_conductivity="%s",\n'
        '  solve_liquid="%s",\n'
        '  solve_grain_growth="%s",\n'
        '  forcing_input_file="%s",\n'
        '  initial_input_file="%s",\n'
        '  time_step=%f,\n'
        '  spinup=%f\n'
        ' /\n'
    ) % (
        solve_density,
        solve_temperature,
        solve_heat_capacity,
        solve_thermal_conductivity,
        solve_liquid,
        solve_grain_growth,
        forcing_input_file,
        initial_input_file,
        time_step,
        spinup,
    )

    f = open('tfm.conf', 'w')
    f.write(conf)
    f.close()



if __name__ == '__main__':
    
    # download data
    urllib.request.urlretrieve(INIT_DENS_URL, INIT_DENS_FILE)
    urllib.request.urlretrieve(INIT_TEMP_URL, INIT_TEMP_FILE)
    urllib.request.urlretrieve(INIT_AGE_URL, INIT_AGE_FILE)
    urllib.request.urlretrieve(FORCING_TEMP_URL, FORCING_TEMP_FILE)
    urllib.request.urlretrieve(FORCING_ACC_URL, FORCING_ACC_FILE)

    # data import
    init_dens = data_import(
        INIT_DENS_FILE,
        cols=INIT_DENS_COLS,
    )[::-1]
    init_temp = data_import(
        INIT_TEMP_FILE,
        cols=INIT_TEMP_COLS,
    )[::-1]
    init_age = data_import(
        INIT_AGE_FILE,
        cols=INIT_AGE_COLS,
    )[::-1]
    forcing_temp = data_import(
        FORCING_TEMP_FILE,
        cols=FORCING_TEMP_COLS,
        dtype=str,
    )[::-1]
    forcing_acc = data_import(
        FORCING_ACC_FILE,
        cols=FORCING_ACC_COLS,
    )[::-1]


    # data conversation initial density
    init_dens[:,0] *= -1.0

    # data conversation initial temperature
    init_temp[:,0] *= -1.0
    init_temp[:,1] += CELSIUS2KELVIN

    # data conversation initial age
    init_age[:,0] *= -1.0
    init_age[:,1] = -(init_age[:,1] - init_age[-1,1]) * SECONDS_YEAR

    # data conversations forcing temperature
    forcing_temp = forcing_temp[1:]

    forcing_temp_date = [date.split('-') for date in forcing_temp[:,0]]
    forcing_temp_date = [
        float(date[0]) + float(date[1]) / 12.0 for date in forcing_temp_date
    ]
    forcing_temp_date = fix_temperature_dates(forcing_temp_date)
    forcing_temp_temp = forcing_temp[:,1].astype(float) + CELSIUS2KELVIN

    forcing_temp = np.array([forcing_temp_date, forcing_temp_temp]).transpose()

    # data conversations forcing accumulation
    forcing_acc[:,1] = forcing_acc[:,1] / 1000.0 / SECONDS_YEAR


    # interpolate inital temperature and age data to inital density data
    min_depth = max([
        min(init_dens[:,0]),
        min(init_temp[:,0]),
        min(init_age[:,0]),
    ])
    max_depth = min([
        max(init_dens[:,0]),
        max(init_temp[:,0]),
        max(init_age[:,0]),
    ])

    mask = np.logical_and(
        (init_dens[:,0] >= min_depth),
        (init_dens[:,0] <= max_depth),
    )
    init_dens = init_dens[mask,:]

    mask = np.logical_and(
        (init_temp[:,0] >= min_depth),
        (init_temp[:,0] <= max_depth),
    )
    init_temp = init_temp[mask,:]
    interp_temp = np.interp(
        init_dens[:,0],
        init_temp[:,0], init_temp[:,1],
    )

    mask = np.logical_and(
        (init_age[:,0] >= min_depth),
        (init_age[:,0] <= max_depth),
    )
    init_age = init_age[mask,:]
    interp_age = np.interp(
        init_dens[:,0],
        init_age[:,0], init_age[:,1],
    )


    # interpolate accumulation data to temperature data
    min_time = max([min(forcing_temp[:,0]), min(forcing_acc[:,0])])
    max_time = min([max(forcing_temp[:,0]), max(forcing_acc[:,0])])

    mask = np.logical_and(
        (forcing_temp[:,0] >= min_time),
        (forcing_temp[:,0] <= max_time),
    )
    forcing_temp = forcing_temp[mask,:]

    mask = np.logical_and(
        (forcing_acc[:,0] >= min_time),
        (forcing_acc[:,0] <= max_time),
    )
    forcing_acc = forcing_acc[mask,:]

    interp_acc = np.interp(
        forcing_temp[:,0],
        forcing_acc[:,0], forcing_acc[:,1],
    )


    # assemble init data
    init_depth        = init_dens[:,0]
    init_dens         = init_dens[:,1]
    init_temp         = interp_temp
    init_grain_radius = np.ones_like(init_depth) * INIT_GRAIN_RADIUS
    init_liquid_water = np.ones_like(init_depth) * INIT_LIQUID_WATER
    init_age          = interp_age

    init = np.array([
        init_depth,
        init_dens,
        init_temp,
        init_grain_radius,
        init_liquid_water,
        init_age,
    ]).transpose()

    np.savetxt(INITIAL_INPUT_FILE, init)


    # assemble final forcing data
    time_step = round(np.mean(abs(np.diff(forcing_temp[:,0]))) * SECONDS_YEAR)

    forcing_world_time   = forcing_temp[:,0]
    forcing_model_time   = np.arange(len(forcing_world_time)) * time_step
    forcing_temp         = forcing_temp[:,1]
    forcing_dens         = np.ones_like(interp_acc) * FORCING_DENSITY
    forcing_solid_acc    = interp_acc
    forcing_liquid_acc   = np.ones_like(interp_acc) * FORCING_LIQUID_ACC
    forcing_grain_radius = np.ones_like(interp_acc) * FORCING_GRAIN_RADIUS

    forcing = np.array([
        forcing_world_time,
        forcing_model_time,
        forcing_temp,
        forcing_dens,
        forcing_solid_acc,
        forcing_liquid_acc,
        forcing_grain_radius,
    ]).transpose()

    np.savetxt(FORCING_INPUT_FILE, forcing)


    # configuration file generation
    generate_config(
        SOLVE_DENSITY,
        SOLVE_TEMPERATURE,
        SOLVE_HEAT_CAPACITY,
        SOLVE_THERMAL_CONDUCTIVITY,
        SOLVE_LIQUID,
        SOLVE_GRAIN_GROWTH,
        FORCING_INPUT_FILE,
        INITIAL_INPUT_FILE,
        time_step,
        SPINUP,
    )
