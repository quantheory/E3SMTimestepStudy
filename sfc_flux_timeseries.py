#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

START_DAY = 1
END_DAY = 1
TIME_STEP = 1800
assert 86400 % TIME_STEP == 0, "cannot fit even number of time steps in day"
times_per_day = 86400 // TIME_STEP

CASE_NAMES = [
    "timestep_ctrl" # SPS: REPLACE THIS
]
SHORT_CASE_NAMES = [
    "CTRL",
]
STYLES = {
    "CTRL": ('k', '-'),
}

OUTPUT_ROOT = "/p/lustre2/santos36/ACME" # SPS: REPLACE THIS

OUTPUT_DIRS = ["{}/{}/run/".format(OUTPUT_ROOT, case)
               for case in CASE_NAMES]

suffix = "_fluxdiag"

log_file = open("sfc_timeseries_log{}.txt".format(suffix), 'w')

out_file_template = "{}.cam.h0.0001-01-{}-{}.nc"

def day_str(day):
    "Given an integer day, return the 2-digit day string used for file names."
    return "{:02d}".format(day)

def time_str(time):
    "Given an integer time in seconds, return the 5-digit string used in file names."
    return "{:05d}".format(time)

def get_out_file_name(icase, day, time):
    """Given a case index, day, and time, return CAM header file name."""
    return join(OUTPUT_DIRS[icase],
                out_file_template.format(CASE_NAMES[icase],
                                         day_str(day), time_str(time)))

first_file_name = get_out_file_name(0, 1, 0)
first_file = nc4.Dataset(first_file_name, 'r')
ncol = len(first_file.dimensions['ncol'])
nlev = len(first_file.dimensions['lev'])
lat = first_file['lat'][:]
lon = first_file['lon'][:]
lev = first_file['lev'][:]

# Find columns in box over South America.
min_lat = -20.
max_lat = 10.
min_lon = 280.
max_lon = 315.
column_set = set()
for i in range(ncol):
    if min_lon <= lon[i] <= max_lon and min_lat <= lat[i] <= max_lat:
        column_set.add(i)

first_file.close()

ncol_sa = len(column_set)

# Plot wind issues over SA.
from mpl_toolkits import basemap
bmap = basemap.Basemap(lon_0=180.)

#TIME_CHECK = 57600
#TIME_CHECK = 79200
TIME_CHECK = 0
#TIME_CHECK = 54000
#TIME_CHECK = 12*3600

DAY_CHECK = 2

CASE_CHECK = 0

LEVEL = nlev - 1

case = SHORT_CASE_NAMES[CASE_CHECK]
time_increment = TIME_STEP

plot_box = [min_lon, max_lon, min_lat, max_lat]
#plot_box = [285., 297., 1., 10.]

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK)
mid_file = nc4.Dataset(mid_file_name, 'r')

u1 = mid_file['U'][0,LEVEL,:]
v1 = mid_file['V'][0,LEVEL,:]

plt.quiver(lon, lat, u1, v1,
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow1{}.png".format(suffix))
plt.close()

mid_file.close()

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK + time_increment)
mid_file = nc4.Dataset(mid_file_name, 'r')

u2 = mid_file['U'][0,LEVEL,:]
v2 = mid_file['V'][0,LEVEL,:]

plt.quiver(lon, lat, u2, v2,
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow2{}.png".format(suffix))
plt.close()

mid_file.close()

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK + 2*time_increment)
mid_file = nc4.Dataset(mid_file_name, 'r')

u3 = mid_file['U'][0,LEVEL,:]
v3 = mid_file['V'][0,LEVEL,:]

plt.quiver(lon, lat, u3, v3,
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow3{}.png".format(suffix))
plt.close()

mid_file.close()

ud2 = u1 - 2*u2 + u3
vd2 = v1 - 2*v2 + v3

# Find column with large oscillations in wind speed
# To force focus on a particular column, just set the number here.
ifocus = -1

if ifocus == -1:
    print("Searching for oscillatory point.", file=log_file, flush=True)
    maxd2 = 0.
    for icol in column_set:
        d2 = np.sqrt(ud2[icol]*ud2[icol] + vd2[icol]*vd2[icol])
        if d2 > maxd2:
            maxd2 = d2
            ifocus = icol

    assert ifocus >= 0, "no focus column found"

print("Worst oscillations at column ", ifocus, " at lat = ",
      lat[ifocus], ", lon = ", lon[ifocus], file=log_file, flush=True)

plt.scatter(lon[ifocus], lat[ifocus])

plt.quiver(lon, lat, ud2, vd2,
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_D2_arrow{}.png".format(suffix))
plt.close()

variables = [
    {'name': 'RELHUM', 'units': r'%', 'ndim': 2},
    {'name': 'CLDLIQ', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'LHFLX', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'SHFLX', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'TS', 'units': r'$K$', 'ndim': 1},
    {'name': 'T', 'units': r'$K$', 'ndim': 2},
    {'name': 'Q', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'U', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'V', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'U10', 'units': r'$m/s$', 'ndim': 1},
    {'name': 'PS', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUX', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUY', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUGWX', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUGWY', 'units': r'$Pa$', 'ndim': 1},
]

derived_variables = [
    {'name': 'TAU', 'units': r'$Pa$', 'ndim': 1,
     'depends': ['TAUX', 'TAUY'],
     'calc': (lambda var_dict: np.sqrt(var_dict['TAUX']**2 + var_dict['TAUY']**2)),
    },
]

# Check that dependencies are satisfied.
var_names = [var['name'] for var in variables]
for derived in derived_variables:
    for depend in derived['depends']:
        assert depend in var_names

ncases = len(CASE_NAMES)
ntimes = (END_DAY - START_DAY + 1) * times_per_day + 1

out_vars = {}
for icase in range(ncases):
    case = SHORT_CASE_NAMES[icase]
    print("Processing case ", case)
    case_times_per_day = times_per_day
    case_ntimes = ntimes
    case_time_step = TIME_STEP
    out_vars[case] = {}
    for var in variables:
        out_vars[case][var['name']] = np.zeros((case_ntimes,))
    ita = 0
    for day in range(START_DAY, END_DAY+1):
        for it in range(case_times_per_day):
            out_file_name = get_out_file_name(icase, day, it*case_time_step)
            out_file = nc4.Dataset(out_file_name, 'r')
            for var in variables:
                varname = var['name']
                ndim = var['ndim']
                if ndim == 1:
                    out_vars[case][varname][ita] = out_file[varname][0,ifocus]
                elif ndim == 2:
                    out_vars[case][varname][ita] = out_file[varname][0,LEVEL,ifocus]
                else:
                    assert False, \
                        "don't know what to do with ndim={}".format(ndim)
            out_file.close()
            ita += 1
    # Last file is 0-th time of the next day.
    out_file_name = get_out_file_name(icase, END_DAY+1, 0)
    out_file = nc4.Dataset(out_file_name, 'r')
    for var in variables:
        varname = var['name']
        ndim = var['ndim']
        if ndim == 1:
            out_vars[case][varname][ita] = out_file[varname][0,ifocus]
        elif ndim == 2:
            out_vars[case][varname][ita] = out_file[varname][0,LEVEL,ifocus]
        else:
            assert False, \
                "don't know what to do with ndim={}".format(ndim)
    out_file.close()
    # Scale variables
    for var in variables:
        if 'scale' in var:
            out_vars[case][var['name']] *= var['scale']
    # Calculate derived variables
    for derived in derived_variables:
        out_vars[case][derived['name']] = derived['calc'](out_vars[case])

# Assumes Venezuelan time.
TIME_OFFSET = 4.
times = np.linspace(0., TIME_STEP*(ntimes - 1) / 3600., ntimes) - TIME_OFFSET
for var in variables + derived_variables:
    name = var['name']
    for icase in range(ncases):
        case = SHORT_CASE_NAMES[icase]
        case_times = times
        plt.plot(case_times, out_vars[case][name], color=STYLES[case][0],
                 linestyle=STYLES[case][1])
    plt.axis('tight')
    plt.xlabel("Time (UTC-4:00)")
    # Bad hard-coding!
    if START_DAY - END_DAY + 1 == 1:
        plt.xticks(np.linspace(-3., 18., 8),
                   ["2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800"])
    elif START_DAY - END_DAY + 1 == 2:
        plt.xticks(np.linspace(-3., 42., 16),
                   ["2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800",
                    "2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800"])
    plt.grid(True)
    if 'display' in var:
        dname = var['display']
    else:
        dname = name
    plt.ylabel("{} ({})".format(dname, var['units']))
    plt.savefig("{}_time{}.png".format(name, suffix))
    plt.close()

log_file.close()
