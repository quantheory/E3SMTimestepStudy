#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str, time_str

CASENAME = "timestep_CLUBB_MG2_10s"

NUM_DAYS = 1
TIME_STEP = 1800
assert 86400 % TIME_STEP == 0, "cannot fit even number of time steps in day"
times_per_day = 86400 // TIME_STEP

CASE_NAMES = [
#    "timestep_ctrl",
#    "timestep_MG2_10s",
#    "timestep_CLUBB_10s_MG2_10s",
#    "timestep_CLUBB_MG2_10s",
#    "timestep_CLUBB_MG2_10s_ftype1",
#    "timestep_all_10s",
#    "timestep_dyn_10s",
    "timestep_presaer_ctrl",
    "timestep_presaer_CLUBB_MG2_10s",
    "timestep_presaer_CLUBB_MG2_10s_ZM_10s",
    "timestep_presaer_cld_10s",
    "timestep_presaer_cld_10s_ftype1",
    "timestep_presaer_all_10s",
]
SHORT_CASE_NAMES = [
#    "CTRL",
#    "MICRO10",
#    "CLUBB10MICRO10",
#    "CLUBBMICRO10",
#    "CLUBBMICRO10FTYPE1",
#    "ALL10",
#    "DYN10",
    "CTRLPA",
    "CLUBBMICRO10PA",
    "CLUBBMICRO10ZM10PA",
    "CLD10PA",
    "CLD10FTYPE1PA",
    "ALL10PA",
]
STYLES = {
    "CTRL": ('k', '-'),
    "MICRO10": ('r', '-'),
    "CLUBB10MICRO10": ('maroon', '-'),
    "CLUBBMICRO10": ('indigo', '-'),
    "CLUBBMICRO10FTYPE1": ('indigo', ':'),
    "ALL10": ('dimgrey', '-'),
    "DYN10": ('y', '-'),
    "CTRLPA": ('k', '-'),
    "CLUBBMICRO10PA": ('indigo', '-'),
    "CLUBBMICRO10ZM10PA": ('saddlebrown', '-'),
    "CLD10PA": ('slateblue', '-'),
    "CLD10FTYPE1PA": ('slateblue', ':'),
    "ALL10PA": ('dimgrey', '-'),
}
OUTPUT_DIRS = ["/p/lustre2/santos36/ACME/{}/run/".format(case)
               for case in CASE_NAMES]

suffix = ""

log_file = open("plot_water_budget_col_log{}.txt".format(suffix), 'w')

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
ilev = first_file['ilev'][:]

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
column_list = sorted(list(column_set))

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

#ifocus = -1

if ifocus == -1:
    # Look at precip in a particular run.
    print("Searching for largest average precipitation.", file=log_file)
    itest = 1
    precl_total = np.zeros((ncol_sa,))
    for day in range(1, NUM_DAYS+1):
        for it in range(times_per_day):
            test_file_name = get_out_file_name(itest, day, it*TIME_STEP)
            test_file = nc4.Dataset(test_file_name, 'r')
            for icol in range(ncol_sa):
                precl_total[icol] += test_file['PRECL'][0,column_list[icol]]
            test_file.close()

    test_file_name = get_out_file_name(itest, NUM_DAYS+1, 0)
    test_file = nc4.Dataset(test_file_name, 'r')
    for icol in range(ncol_sa):
        precl_total[icol] += test_file['PRECL'][0,column_list[icol]]
    test_file.close()

    precl_max = 0.
    for icol in range(ncol_sa):
        if precl_max < precl_total[icol]:
            precl_max = precl_total[icol]
            ifocus = column_list[icol]

    assert ifocus != -1, "Cloud liquid max difference not found!"

print("Difference maximized at column ", ifocus, " at lat = ",
      lat[ifocus], ", lon = ", lon[ifocus], file=log_file)

first_file_name = get_out_file_name(0, 1, 0)
first_file = nc4.Dataset(first_file_name, 'r')
p0 = first_file['P0']
ps = first_file['PS'][0,ifocus]
hyam = first_file['hyam'][:]
hybm = first_file['hybm'][:]
p = hyam * p0 + hybm * ps
first_file.close()

variables = [
    {'name': 'RELHUM', 'units': r'%', 'ndim': 2},
    {'name': 'CLDLIQ', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'CLDICE', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'RAINQM', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'CMELIQ', 'units': r'$g/kg/d$', 'ndim': 2, 'scale': 86.4e6},
    {'name': 'PRAO', 'units': r'$g/kg/d$', 'ndim': 2, 'scale': 86.4e6},
    {'name': 'PRCO', 'units': r'$g/kg/d$', 'ndim': 2, 'scale': 86.4e6},
    {'name': 'QCSEDTEN', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'QRSEDTEN', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'EVAPPREC', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'T', 'units': r'$K$', 'ndim': 2},
    {'name': 'Q', 'units': r'$g/kg$', 'ndim': 2, 'scale': 1000.},
    {'name': 'U', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'V', 'units': r'$m/s$', 'ndim': 2},
    {'name': 'OMEGA', 'units': r'$Pa/s$', 'ndim': 2},
    {'name': 'CLOUD', 'units': r'fraction', 'ndim': 2},
    {'name': 'DPDLFLIQ', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'DPDLFICE', 'units': r'$kg/kg/s$', 'ndim': 2},
    {'name': 'QRL', 'units': r'$K/s$', 'ndim': 2},
    {'name': 'QRS', 'units': r'$K/s$', 'ndim': 2},
    {'name': 'Z3', 'units': r'$m$', 'ndim': 2},
    {'name': 'QCSEVAP', 'units': r'$g/kg/d$', 'ndim': 2, 'scale': 86.4e6},
    {'name': 'QISEVAP', 'units': r'$g/kg/d$', 'ndim': 2, 'scale': 86.4e6},
]

def calc_rho(t):
    rho = np.zeros(t.shape)
    ntimes = t.shape[1]
    for i in range(ntimes):
        rho[:,i] = p / (287.058 * t[:,i])
    return rho

derived_variables = [
    {'name': 'LWC', 'units': r'$g/m^3$', 'ndim': 2,
     'depends': ['CLDLIQ', 'T'],
     'calc': (lambda var_dict: var_dict['CLDLIQ'] * calc_rho(var_dict['T'])),
    },
    {'name': 'IWC', 'units': r'$g/m^3$', 'ndim': 2,
     'depends': ['CLDICE', 'T'],
     'calc': (lambda var_dict: var_dict['CLDICE'] * calc_rho(var_dict['T'])),
    },
    {'name': 'RAINPROD', 'units': r'$g/kg/d$', 'ndim': 2,
     'depends': ['PRCO', 'PRAO'],
     'calc': (lambda var_dict: var_dict['PRCO'] + var_dict['PRAO']),
    },
]

# Check that dependencies are satisfied.
var_names = [var['name'] for var in variables]
for derived in derived_variables:
    for depend in derived['depends']:
        assert depend in var_names

ncases = len(CASE_NAMES)
ntimes = NUM_DAYS * times_per_day + 1

out_vars = {}
for icase in range(ncases):
    case = SHORT_CASE_NAMES[icase]
    print("Processing case ", case)
    out_vars[case] = {}
    for var in variables:
        out_vars[case][var['name']] = np.zeros((nlev, ntimes))
    ita = 0
    for day in range(1, NUM_DAYS+1):
        for it in range(times_per_day):
            out_file_name = get_out_file_name(icase, day, it*TIME_STEP)
            out_file = nc4.Dataset(out_file_name, 'r')
            for var in variables:
                varname = var['name']
                ndim = var['ndim']
                if ndim == 2:
                    out_vars[case][varname][:,ita] = out_file[varname][0,:,ifocus]
                else:
                    assert False, \
                        "don't know what to do with ndim={}".format(ndim)
            out_file.close()
            ita += 1
    # Last file is 0-th time of the next day.
    out_file_name = get_out_file_name(icase, NUM_DAYS+1, 0)
    out_file = nc4.Dataset(out_file_name, 'r')
    for var in variables:
        varname = var['name']
        ndim = var['ndim']
        if ndim == 2:
            out_vars[case][varname][:,ita] = out_file[varname][0,:,ifocus]
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

PLOT_TOP = 100.
itop = 0
for level in ilev:
    if level > PLOT_TOP:
        break
    itop += 1

plot_ilev = ilev[itop:]

# Assumes Venezuelan time.
TIME_OFFSET = 4.
times = np.linspace(0., TIME_STEP*(ntimes - 1) / 3600., ntimes) - TIME_OFFSET
for var in variables + derived_variables:
    name = var['name']
    clim_val = [1.e36, -1.e36]
    for icase in range(ncases):
        clim_val[0] = min(clim_val[0], out_vars[case][name][itop:,1:].min())
        clim_val[1] = max(clim_val[1], out_vars[case][name][itop:,1:].max())

    for icase in range(ncases):
        case = SHORT_CASE_NAMES[icase]
        plt.pcolor(times, plot_ilev, out_vars[case][name][itop:,1:])
        plt.axis('tight')
        plt.xlabel("Time (hr)")
        # Bad hard-coding!
        if NUM_DAYS == 1:
            plt.xticks(np.linspace(-3., 18., 8),
                       ["2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800"])
        elif NUM_DAYS == 2:
            plt.xticks(np.linspace(-3., 42., 16),
                       ["2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800",
                        "2100", "0000", "0300", "0600", "0900", "1200", "1500", "1800"])
        plt.grid(True, axis='x')
        ax = plt.gca()
        ylim = ax.get_ylim()
        ax.set_ylim([ylim[1], ylim[0]])
        plt.ylabel("Pressure (hPa)")
        plt.colorbar()
        plt.clim(clim_val[0], clim_val[1])
        plt.savefig("{}_{}_time_col{}.png".format(name, case, suffix))
        plt.close()

htime = (15 * 3600) // TIME_STEP

for icase in range(ncases):
    case = SHORT_CASE_NAMES[icase]
    hodo_winds = np.zeros((2, 11))
    next_km = 1
    base_z3 = out_vars[case]['Z3'][nlev-1,htime]
    for jlev in range(nlev-1, -1, -1):
        this_z3 = out_vars[case]['Z3'][jlev,htime]
        next_z3 = out_vars[case]['Z3'][jlev-1,htime]
        if next_z3 > (1000.*next_km + base_z3):
            weight = (1000.*next_km + base_z3 - this_z3) / (next_z3 - this_z3)
            u = out_vars[case]['U'][jlev-1,htime]*weight + \
                out_vars[case]['U'][jlev,htime]*(1.-weight)
            v = out_vars[case]['V'][jlev-1,htime]*weight + \
                out_vars[case]['V'][jlev,htime]*(1.-weight)
            hodo_winds[0,next_km] = u
            hodo_winds[1,next_km] = v
            next_km += 1
        if next_km == 11:
            break
    plt.plot(hodo_winds[0,:], hodo_winds[1,:])
    plt.savefig("hodo_{}{}.png".format(case, suffix))
    plt.close()

log_file.close()
