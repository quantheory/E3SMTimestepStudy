#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str, time_str

START_DAY = 1
END_DAY = 1
TIME_STEP = 1800
assert 86400 % TIME_STEP == 0, "cannot fit even number of time steps in day"
times_per_day = 86400 // TIME_STEP

CASE_NAMES = [
    "timestep_ctrl",
#    "timestep_MG2_10s",
    "timestep_CLUBB_10s_MG2_10s",
    "timestep_CLUBB_MG2_10s",
#    "timestep_CLUBB_MG2_10s_ftype1",
    "timestep_all_10s",
#    "timestep_all_300s",
#    "timestep_all_900s",
#    "timestep_ctrl_notms",
#    "timestep_ctrl_niter_change",
#    "timestep_ctrl_niter_change_smooth",
#    "timestep_smooth14",
#    "timestep_smooth35",
#    "timestep_300s_niter_change",
#    "timestep_dyn_10s",
#    "timestep_presaer_ctrl",
#    "timestep_presaer_CLUBB_MG2_10s",
#    "timestep_presaer_CLUBB_MG2_10s_ZM_10s",
#    "timestep_presaer_cld_10s",
#    "timestep_presaer_cld_10s_ftype1",
#    "timestep_presaer_all_10s",
#    "timestep_presaer_all_10s_lower_tau",
#    "timestep_presaer_all_10s_cld_5s",
]
SHORT_CASE_NAMES = [
    "CTRL",
#    "MICRO10",
    "CLUBB10MICRO10",
    "CLUBBMICRO10",
#    "CLUBBMICRO10FTYPE1",
    "ALL10",
#    "ALL300",
#    "ALL900",
#    "CTRLFLXAVG",
#    "CTRLNITERS",
#    "CTRLNITERSSMOOTH35",
#    "SMOOTH14",
#    "SMOOTH35",
#    "ALL300NITERS",
#    "DYN10",
#    "CTRLPA",
#    "CLUBBMICRO10PA",
#    "CLUBBMICRO10ZM10PA",
#    "CLD10PA",
#    "CLD10FTYPE1PA",
#    "ALL10PA",
#    "ALL10PALT",
#    "ALL10CLD5PA",
]
STYLES = {
    "CTRL": ('k', '-'),
    "MICRO10": ('r', '-'),
    "CLUBB10MICRO10": ('maroon', '-'),
    "CLUBBMICRO10": ('indigo', '-'),
    "CLUBBMICRO10FTYPE1": ('indigo', ':'),
    "ALL10": ('dimgrey', '-'),
#    "ALL300": ('dimgrey', ':'),
    "ALL300": ('g', '-'),
    "ALL900": ('k', '--'),
    "CTRLFLXAVG": ('k', ':'),
    "CTRLNITERS": ('orange', '-'),
    "CTRLNITERSSMOOTH35": ('orange', '--'),
    "SMOOTH14": ('r', '-'),
    "SMOOTH35": ('r', '-'),
    "ALL300NITERS": ('b', '-'),
    "DYN10": ('y', '-'),
    "CTRLPA": ('k', '-'),
    "CLUBBMICRO10PA": ('indigo', '-'),
    "CLUBBMICRO10ZM10PA": ('saddlebrown', '-'),
    "CLD10PA": ('slateblue', '-'),
    "CLD10FTYPE1PA": ('slateblue', ':'),
    "ALL10PA": ('dimgrey', '-'),
    "ALL10PALT": ('dimgrey', '-.'),
    "ALL10CLD5PA": ('slateblue', '-.'),
}
OUTPUT_DIRS = ["/p/lustre2/santos36/ACME/{}/run/".format(case)
               for case in CASE_NAMES]

suffix = ""

log_file = open("plot_water_budget_box_log{}.txt".format(suffix), 'w')

out_file_template = "{}.cam.h0.0001-01-{}-{}.nc"

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

# Max diff in CLDLIQ at surface for CLUBBMICRO10
#ifocus = 28970
# Max CLDLIQ at surface for CLUBBMICRO10
#ifocus = 27898
# Max precipitation in CTRL
#ifocus = 29215
# Max precipitation in CLUBBMICRO10 and CLUBBMICRO10PA
ifocus = 29488
# Max precipitation in ALL10
#ifocus = 29227
# Large oscillations at 1800
#ifocus = 29520

#ifocus = -1

if ifocus == -1:
    # Look at lowest level CLDLIQ, to see where fog difference is most intense
    # between two cases.
    print("Searching for largest fog difference.", file=log_file, flush=True)
    ictrl = 0
    itest = 2
    cldliq_maxdiff = 0.
    ifocus = -1
    for day in range(START_DAY, END_DAY+1):
        for it in range(times_per_day):
            ctrl_file_name = get_out_file_name(ictrl, day, it*TIME_STEP)
            ctrl_file = nc4.Dataset(ctrl_file_name, 'r')
            test_file_name = get_out_file_name(itest, day, it*TIME_STEP)
            test_file = nc4.Dataset(test_file_name, 'r')
            for icol in column_set:
                cldliq_diff = test_file['CLDLIQ'][0,nlev-1,icol] - \
                              ctrl_file['CLDLIQ'][0,nlev-1,icol]
                if cldliq_diff > cldliq_maxdiff:
                    cldliq_maxdiff = cldliq_diff
                    ifocus = icol
            ctrl_file.close()
            test_file.close()

    ctrl_file_name = get_out_file_name(ictrl, END_DAY+1, 0)
    ctrl_file = nc4.Dataset(ctrl_file_name, 'r')
    test_file_name = get_out_file_name(itest, END_DAY+1, 0)
    test_file = nc4.Dataset(test_file_name, 'r')
    for icol in column_set:
        cldliq_diff = test_file['CLDLIQ'][0,nlev-1,icol] - \
                      ctrl_file['CLDLIQ'][0,nlev-1,icol]
        if cldliq_diff > cldliq_maxdiff:
            cldliq_maxdiff = cldliq_diff
            ifocus = icol
    ctrl_file.close()
    test_file.close()

    assert ifocus != -1, "Cloud liquid max difference not found!"

#print("Difference maximized at column ", ifocus, " at lat = ",
#      lat[ifocus], ", lon = ", lon[ifocus], file=log_file, flush=True)

variables = [
    {'name': 'RELHUM', 'units': r'%', 'ndim': 2},
    {'name': 'CLDLIQ', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'RAINQM', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'LHFLX', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'SHFLX', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'CMELIQ', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'PRAO', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'PRCO', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'QCSEDTEN', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'QRSEDTEN', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'PRECL', 'units': r'$mm/day$', 'ndim': 1, 'scale': 1000.*86400.},
    {'name': 'PRECC', 'units': r'$mm/day$', 'ndim': 1, 'scale': 1000.*86400.},
    {'name': 'EVAPPREC', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'T', 'units': r'$K$', 'ndim': 2},
    {'name': 'Q', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'U', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'V', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'OMEGA', 'units': r'$Pa/s$', 'ndim': 2},
    {'name': 'OMEGA500', 'units': r'$Pa/s$', 'ndim': 1, 'display': r'$\omega_{500}$'},
    {'name': 'U10', 'units': r'$m/s$', 'ndim': 1},
    {'name': 'PSL', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'FSDS', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'SWCF', 'units': r'$W/m^2$', 'ndim': 1},
    {'name': 'CLDLOW', 'units': r'fraction', 'ndim': 1},
    {'name': 'CAPE', 'units': r'$J/kg$', 'ndim': 1},
    {'name': 'TAUX', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUY', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUGWX', 'units': r'$Pa$', 'ndim': 1},
    {'name': 'TAUGWY', 'units': r'$Pa$', 'ndim': 1},
]

derived_variables = [
    {'name': 'PRECT', 'units': r'$mm/day$', 'ndim': 1,
     'depends': ['PRECL', 'PRECC'],
     'calc': (lambda var_dict: var_dict['PRECL'] + var_dict['PRECC']),
     'display': "Total Precipitation",
    },
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

LEVEL = nlev - 1

out_vars = {}
for icase in range(ncases):
    case = SHORT_CASE_NAMES[icase]
    print("Processing case ", case)
    if case == "ALL900":
        case_times_per_day = times_per_day * 2
        case_ntimes = ntimes*2 - 1
        case_time_step = TIME_STEP // 2
    elif case == "ALL300" or case == "ALL300NITERS":
        case_times_per_day = times_per_day * 6
        case_ntimes = ntimes*6 - 5
        case_time_step = TIME_STEP // 6
    else:
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
times900 = np.linspace(0., TIME_STEP*(ntimes - 1) / 3600., 2*ntimes - 1) - TIME_OFFSET
times300 = np.linspace(0., TIME_STEP*(ntimes - 1) / 3600., 6*ntimes - 5) - TIME_OFFSET
print_vars = {"PRECT", "PRECL", "PRECC", "OMEGA500", "CAPE"}
print("time", ",", ",".join(str(time % 24) for time in times), sep="",
      file=log_file, flush=True)
for var in variables + derived_variables:
    name = var['name']
    for icase in range(ncases):
        case = SHORT_CASE_NAMES[icase]
        if case == "ALL900":
            case_times = times900
        elif case == "ALL300" or case == "ALL300NITERS":
            case_times = times300
        else:
            case_times = times
        if name in print_vars:
            print(name, ",", case, ",", ",".join(str(x) for x in out_vars[case][name]),
                  sep="", file=log_file, flush=True)
        plt.plot(case_times, out_vars[case][name], color=STYLES[case][0],
                 linestyle=STYLES[case][1])
    plt.axis('tight')
    plt.xlabel("Time (hr)")
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
    plt.savefig("{}_time_box{}.png".format(name, suffix))
    plt.close()

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

case = SHORT_CASE_NAMES[CASE_CHECK]
if case == "ALL900":
    time_increment = TIME_STEP // 2
elif case == "ALL300" or case == "ALL300NITERS":
    time_increment = TIME_STEP // 6
else:
    time_increment = TIME_STEP

#plot_box = [min_lon, max_lon, min_lat, max_lat]
plot_box = [285., 297., 1., 10.]

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK)
mid_file = nc4.Dataset(mid_file_name, 'r')

u1 = mid_file['U'][0,LEVEL,:]
v1 = mid_file['V'][0,LEVEL,:]

#plt.scatter(lon, lat, c=mid_file['PSL'][0,:])
#plt.colorbar()
plt.quiver(lon, lat, mid_file['U'][0,LEVEL,:], mid_file['V'][0,LEVEL,:],
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow_box1{}.png".format(suffix))
plt.close()

plt.scatter(lon, lat, c=mid_file['OMEGA500'][0,:])
bmap.drawcoastlines()
plt.axis(plot_box)
plt.colorbar()
plt.savefig("OMEGA500_box1{}.png".format(suffix))
plt.close()

mid_file.close()

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK + time_increment)
mid_file = nc4.Dataset(mid_file_name, 'r')

u2 = mid_file['U'][0,LEVEL,:]
v2 = mid_file['V'][0,LEVEL,:]

plt.quiver(lon, lat, mid_file['U'][0,LEVEL,:], mid_file['V'][0,LEVEL,:],
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow_box2{}.png".format(suffix))
plt.close()

plt.scatter(lon, lat, c=mid_file['OMEGA500'][0,:])
bmap.drawcoastlines()
plt.axis(plot_box)
plt.colorbar()
plt.savefig("OMEGA500_box2{}.png".format(suffix))
plt.close()

mid_file.close()

mid_file_name = get_out_file_name(CASE_CHECK, DAY_CHECK, TIME_CHECK + 2*time_increment)
mid_file = nc4.Dataset(mid_file_name, 'r')

u3 = mid_file['U'][0,LEVEL,:]
v3 = mid_file['V'][0,LEVEL,:]

plt.quiver(lon, lat, mid_file['U'][0,LEVEL,:], mid_file['V'][0,LEVEL,:],
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_arrow_box3{}.png".format(suffix))
plt.close()

plt.scatter(lon, lat, c=mid_file['OMEGA500'][0,:])
bmap.drawcoastlines()
plt.axis(plot_box)
plt.colorbar()
plt.savefig("OMEGA500_box3{}.png".format(suffix))
plt.close()

mid_file.close()

plt.scatter(lon[ifocus], lat[ifocus])

plt.quiver(lon, lat, u1 - 2*u2 + u3, v1 - 2*v2 + v3,
           scale=100., scale_units='height', angles='xy')
bmap.drawcoastlines()
plt.axis(plot_box)
plt.savefig("UV_D2_arrow_box{}.png".format(suffix))
plt.close()

log_file.close()
