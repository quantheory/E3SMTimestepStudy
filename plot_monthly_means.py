#!/usr/bin/env python3

from functools import partial

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

START_YEAR = 1
START_MONTH = 1
END_YEAR = 4
END_MONTH = 2

START_AVG_YEAR = 1
START_AVG_MONTH = 3
END_AVG_YEAR = 4
END_AVG_MONTH = 2

MONTHLY_FILE_LOC="/p/lscratchh/santos36/timestep_monthly_avgs/"

USE_PRESAER=False
TROPICS_ONLY=False
MIDLATITUDES_ONLY=False

assert not (TROPICS_ONLY and MIDLATITUDES_ONLY), \
    "can't do only tropics and only midlatitudes"

nmonths = (END_YEAR - START_YEAR) * 12 - (START_MONTH - 1) + END_MONTH
imonths = list(range(nmonths))
curr_month = START_MONTH
curr_year = START_YEAR
months = []
years = []
for i in range(nmonths):
    months.append(curr_month)
    years.append(curr_year)
    curr_month += 1
    if curr_month > 12:
        curr_month = 1
        curr_year += 1
navgmonths = (END_AVG_YEAR - START_AVG_YEAR) * 12 \
    - (START_AVG_MONTH - 1) + END_AVG_MONTH

suffix = '_y{}m{}-y{}m{}'.format(day_str(START_YEAR),
                                 day_str(START_MONTH),
                                 day_str(END_YEAR),
                                 day_str(END_MONTH))
if USE_PRESAER:
    suffix += '_presaer'
if TROPICS_ONLY:
    suffix += '_tropics'
if MIDLATITUDES_ONLY:
    suffix += '_midlats'

log_file = open("plot_monthly_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", MONTHLY_FILE_LOC, START_MONTH, END_MONTH)
    TEST_CASES = [
    ]
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", MONTHLY_FILE_LOC, START_MONTH, END_MONTH)
    TEST_CASES = [
        E3SMCaseOutput("timestep_all_10s", "ALL10", MONTHLY_FILE_LOC, START_MONTH, END_MONTH),
    ]

case_num = len(TEST_CASES)

rfile0 = nc4.Dataset(REF_CASE.get_monthly_file_name(START_MONTH, START_YEAR), 'r')
nlev = len(rfile0.dimensions['lev'])
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
# For tropics_only cases, just use a weight of 0 for all other cases.
if TROPICS_ONLY:
    lat = rfile0['lat'][:]
    for i in range(ncol):
        if np.abs(lat[i]) > 30.:
            area[i] = 0.
# Same for midlatitudes.
elif MIDLATITUDES_ONLY:
    lat = rfile0['lat'][:]
    for i in range(ncol):
        if np.abs(lat[i]) < 30. or np.abs(lat[i]) > 60.:
            area[i] = 0.
area_sum = area.sum()
weights = area/area_sum
rfile0.close()

def calc_var_stats(ref_case, test_cases, month, year, varnames):
    varnames_read = [name for name in varnames if name != "PRECT" and name != "TAU"]
    if "PRECT" in varnames:
        if "PRECL" not in varnames:
            varnames_read.append("PRECL")
        if "PRECC" not in varnames:
            varnames_read.append("PRECC")
    if "TAU" in varnames:
        if "TAUX" not in varnames:
            varnames_read.append("TAUX")
        if "TAUY" not in varnames:
            varnames_read.append("TAUY")
    ref_time_avg, test_time_avgs, diff_time_avgs = ref_case.compare_monthly_averages(test_cases, month, year, varnames_read)
    if "PRECT" in varnames:
        ref_time_avg["PRECT"] = ref_time_avg["PRECL"] + ref_time_avg["PRECC"]
        for icase in range(case_num):
            test_time_avgs[icase]["PRECT"] = test_time_avgs[icase]["PRECL"] + test_time_avgs[icase]["PRECC"]
            diff_time_avgs[icase]["PRECT"] = diff_time_avgs[icase]["PRECL"] + diff_time_avgs[icase]["PRECC"]
    if "TAU" in varnames:
        ref_time_avg["TAU"] = np.sqrt(ref_time_avg["TAUX"]**2 + ref_time_avg["TAUY"]**2)
        for icase in range(case_num):
            test_time_avgs[icase]["TAU"] = np.sqrt(test_time_avgs[icase]["TAUX"]**2 + test_time_avgs[icase]["TAUY"]**2)
            diff_time_avgs[icase]["TAU"] = test_time_avgs[icase]["TAU"] - ref_time_avg["TAU"]
    ref_avg = dict()
    test_avgs = dict()
    diff_avgs = dict()
    rmses = dict()
    for varname in varnames:
        if varname in vars_3D:
            ref_avg[varname] = np.zeros((nlev,))
            for jlev in range(nlev):
                ref_avg[varname][jlev] = (ref_time_avg[varname][jlev,:] * weights).sum()
        else:
            ref_avg[varname] = (ref_time_avg[varname] * weights).sum()
        test_avgs[varname] = []
        diff_avgs[varname] = []
        rmses[varname] = []
        for i in range(len(test_cases)):
            if varname in vars_3D:
                test_avgs[varname].append(np.zeros((nlev,)))
                diff_avgs[varname].append(np.zeros((nlev,)))
                rmses[varname].append(np.zeros((nlev,)))
                for jlev in range(nlev):
                    test_avgs[varname][-1][jlev] = (test_time_avgs[i][varname][jlev,:] * weights).sum()
                    diff_avgs[varname][-1][jlev] = (diff_time_avgs[i][varname][jlev,:] * weights).sum()
                    rmses[varname][-1][jlev] = np.sqrt((diff_time_avgs[i][varname][jlev,:]**2 * weights).sum())
            else:
                test_avgs[varname].append((test_time_avgs[i][varname] * weights).sum())
                diff_avgs[varname].append((diff_time_avgs[i][varname] * weights).sum())
                rmses[varname].append(np.sqrt((diff_time_avgs[i][varname]**2 * weights).sum()))
            assert np.isclose(diff_avgs[varname][i], test_avgs[varname][i] - ref_avg[varname]).all(), \
                "Problem with diff of variable {} from case {}".format(varname, TEST_CASES[i].short_name)
    return (ref_avg, test_avgs, diff_avgs, rmses)

# Possible ways to extract a 2D section start here:
def identity(x):
    return x

def slice_at(level, x):
    return x[:,level]

def plot_vars_over_time(names, units, scales, log_plot_names):
    ref_means = dict()
    test_means = dict()
    diff_means = dict()
    rmses = dict()
    for name in names:
        if name in vars_3D:
            ref_means[name] = np.zeros((nmonths, nlev))
            test_means[name] = np.zeros((case_num, nmonths, nlev))
            diff_means[name] = np.zeros((case_num, nmonths, nlev))
            rmses[name] = np.zeros((case_num, nmonths, nlev))
        else:
            ref_means[name] = np.zeros((nmonths,))
            test_means[name] = np.zeros((case_num, nmonths))
            diff_means[name] = np.zeros((case_num, nmonths))
            rmses[name] = np.zeros((case_num, nmonths))

    for imonth in range(nmonths):
        month = months[imonth]
        year = years[imonth]
        print("On month: ", month, ", year: ", year, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_var_stats(REF_CASE, TEST_CASES, month, year, names)
        for name in names:
            ref_means[name][imonth] = ref_mean[name]*scales[name]
            for i in range(case_num):
                test_means[name][i,imonth] = test_case_means[name][i]*scales[name]
                diff_means[name][i,imonth] = diff_case_means[name][i]*scales[name]
                rmses[name][i,imonth] = case_rmses[name][i]*scales[name]

    for name in names:
        plot_name = name
        if name in plot_names:
            plot_name = plot_names[name]

        get_2D = identity
        if name == "RELHUM" or name == "Q" or name == "T":
            get_2D = partial(slice_at, nlev-1)

        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        ref_plot_var = get_2D(ref_means[name])
        plot_var(imonths, ref_plot_var, label=REF_CASE.short_name)
        for i in range(case_num):
            test_plot_var = get_2D(test_means[name][i])
            plot_var(imonths,
                     test_plot_var,
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            diff_plot_var = get_2D(diff_means[name][i])
            plot_var(imonths,
                     diff_plot_var,
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("Mean {} difference ({})".format(plot_name, units[name]))
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            rmse_plot_var = get_2D(rmses[name][i])
            plot_var(imonths,
                     rmse_plot_var,
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("{} RMSE ({})".format(plot_name, units[name]))
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        month_weights = [month_days[months[imonth] - 1] for imonth in imonths]
        weight_norm = sum(month_weights)
        month_weights = [weight / weight_norm for weight in month_weights]
        print(name, " has reference mean: ", sum([ref_plot_var[imonth] * month_weights[imonth] for imonth in imonths]),
              file=log_file)
        for i in range(case_num):
            test_plot_var = get_2D(test_means[name][i])
            diff_plot_var = get_2D(diff_means[name][i])
            print(name, " has case ", TEST_CASES[i].short_name, " mean: ", sum([test_plot_var[imonth] * month_weights[imonth] for imonth in imonths]),
                  file=log_file)
            print(name, " has difference mean: ", sum([diff_plot_var[imonth] * month_weights[imonth] for imonth in imonths]),
                  file=log_file)

plot_names = {
    'LWCF': "long wave cloud forcing",
    'SWCF': "short wave cloud forcing",
    'PRECC': "convective precipitation",
    'PRECL': "large scale precipitation",
    'PRECT': "total precipitation",
    'TGCLDIWP': "ice water path",
    'TGCLDLWP': "liquid water path",
    'CLDTOT': "cloud area fraction",
    'CLDLOW': "low cloud area fraction",
    'CLDMED': "medium cloud area fraction",
    'CLDHGH': "high cloud area fraction",
    'LHFLX': "latent heat flux",
    'SHFLX': "sensible heat flux",
    'TAU': "surface wind stress",
    'TS': "surface temperature",
    'PSL': "sea level pressure",
    'OMEGA500': "vertical velocity at 500 mb",
    'U10': "10 meter wind speed",
    'RELHUM': "surface relative humidity",
    'Q': "surface specific humidity",
    'TMQ': "precipitable water",
    'T': "lowest level temperature",
}

units = {
    'LWCF': r'$W/m^2$',
    'SWCF': r'$W/m^2$',
    'PRECC': r'$mm/day$',
    'PRECL': r'$mm/day$',
    'PRECT': r'$mm/day$',
    'TGCLDIWP': r'$kg/m^2$',
    'TGCLDLWP': r'$kg/m^2$',
    'AODABS': r'units?',
    'AODUV': r'units?',
    'AODVIS': r'units?',
    'FLDS': r'$W/m^2$',
    'FLNS': r'$W/m^2$',
    'FLNSC': r'$W/m^2$',
    'FLNT': r'$W/m^2$',
    'FLNTC': r'$W/m^2$',
    'FLUT': r'$W/m^2$',
    'FLUTC': r'$W/m^2$',
    'FSDS': r'$W/m^2$',
    'FSDSC': r'$W/m^2$',
    'FSNS': r'$W/m^2$',
    'FSNSC': r'$W/m^2$',
    'FSNT': r'$W/m^2$',
    'FSNTC': r'$W/m^2$',
    'FSNTOA': r'$W/m^2$',
    'FSNTOAC': r'$W/m^2$',
    'FSUTOA': r'$W/m^2$',
    'FSUTOAC': r'$W/m^2$',
    'CLDTOT': r'fraction',
    'CLDLOW': r'fraction',
    'CLDMED': r'fraction',
    'CLDHGH': r'fraction',
    'OMEGA500': r'Pa/s',
    'LHFLX': r'$W/m^2$',
    'SHFLX': r'$W/m^2$',
    'TAU': r'$N/m^2$',
    'TAUX': r'$N/m^2$',
    'TAUY': r'$N/m^2$',
    'TS': r'$K$',
    'PSL': r'$Pa$',
    'U10': r'$m/s$',
    'RELHUM': r'%',
    'Q': r'$g/kg$',
    'TMQ': r'$kg/m^2$',
    'T': r'$K$',
}
names = list(units.keys())
scales = dict()
for name in names:
    scales[name] = 1.
scales['SWCF'] = -1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.
scales['PRECT'] = 1000.*86400.
scales['Q'] = 1000.

vars_3D = [
    'RELHUM',
    'Q',
    'T',
]

log_plot_names = []#'AODABS', 'AODVIS', 'AODUV']

plot_vars_over_time(names, units, scales, log_plot_names)

log_file.close()
