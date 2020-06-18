#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

START_MONTH = 1
END_MONTH = 12

START_AVG_MONTH = 1
END_AVG_MONTH = 12

MONTHLY_FILE_LOC="/p/lscratchh/santos36/timestep_monthly_avgs/"

USE_PRESAER=False
TROPICS_ONLY=True

months = list(range(START_MONTH, END_MONTH+1))
nmonths = len(months)
navgmonths = END_AVG_MONTH - START_AVG_MONTH + 1

suffix = '_m{}-{}'.format(day_str(START_MONTH), day_str(END_MONTH))
if USE_PRESAER:
    suffix += '_presaer'
if TROPICS_ONLY:
    suffix += '_tropics'

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

rfile0 = nc4.Dataset(REF_CASE.get_monthly_file_name(START_MONTH), 'r')
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
# For tropics_only cases, just use a weight of 0 for all other cases.
if TROPICS_ONLY:
    lat = rfile0['lat'][:]
    for i in range(ncol):
        if np.abs(lat[i]) > 30.:
            area[i] = 0.
area_sum = area.sum()
weights = area/area_sum
rfile0.close()

def calc_2D_var_stats(ref_case, test_cases, month, varnames):
    varnames_read = [name for name in varnames if name != "PRECT"]
    if "PRECT" in varnames:
        if "PRECL" not in varnames:
            varnames_read.append("PRECL")
        if "PRECC" not in varnames:
            varnames_read.append("PRECC")
    ref_time_avg, test_time_avgs, diff_time_avgs = ref_case.compare_monthly_averages(test_cases, month, varnames_read)
    if "PRECT" in varnames:
        ref_time_avg["PRECT"] = ref_time_avg["PRECL"] + ref_time_avg["PRECC"]
        for icase in range(case_num):
            test_time_avgs[icase]["PRECT"] = test_time_avgs[icase]["PRECL"] + test_time_avgs[icase]["PRECC"]
            diff_time_avgs[icase]["PRECT"] = diff_time_avgs[icase]["PRECL"] + diff_time_avgs[icase]["PRECC"]
    ref_avg = dict()
    test_avgs = dict()
    diff_avgs = dict()
    rmses = dict()
    for varname in varnames:
        ref_avg[varname] = (ref_time_avg[varname] * weights).sum()
        test_avgs[varname] = []
        diff_avgs[varname] = []
        rmses[varname] = []
        for i in range(len(test_cases)):
            test_avgs[varname].append((test_time_avgs[i][varname] * weights).sum())
            diff_avgs[varname].append((diff_time_avgs[i][varname] * weights).sum())
            assert np.isclose(diff_avgs[varname][i], test_avgs[varname][i] - ref_avg[varname]), \
                "Problem with diff of variable {} from case {}".format(varname, TEST_CASE_NAMES[i])
            rmses[varname].append(np.sqrt((diff_time_avgs[i][varname]**2 * weights).sum()))
    return (ref_avg, test_avgs, diff_avgs, rmses)

def plot_vars_over_time(names, units, scales, log_plot_names):
    ref_means = dict()
    test_means = dict()
    diff_means = dict()
    rmses = dict()
    for name in names:
        ref_means[name] = np.zeros((nmonths,))
        test_means[name] = np.zeros((case_num, nmonths))
        diff_means[name] = np.zeros((case_num, nmonths))
        rmses[name] = np.zeros((case_num, nmonths))

    for imonth in range(nmonths):
        month = months[imonth]
        print("On month: ", month, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_2D_var_stats(REF_CASE, TEST_CASES, month, names)
        for name in names:
            ref_means[name][imonth] = ref_mean[name]*scales[name]
            for i in range(case_num):
                test_means[name][i,imonth] = test_case_means[name][i]*scales[name]
                diff_means[name][i,imonth] = diff_case_means[name][i]*scales[name]
                rmses[name][i,imonth] = case_rmses[name][i]*scales[name]

    for name in names:
        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        plot_var(months, ref_means[name], label=REF_CASE.short_name)
        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_MONTH
            end_ind = TEST_CASES[i].end_day - START_MONTH + 1
            plot_var(months[start_ind:end_ind],
                     test_means[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("Mean {} ({})".format(name, units[name]))
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_MONTH
            end_ind = TEST_CASES[i].end_day - START_MONTH + 1
            plot_var(months[start_ind:end_ind],
                     diff_means[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("Mean {} difference ({})".format(name, units[name]))
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_MONTH
            end_ind = TEST_CASES[i].end_day - START_MONTH + 1
            plot_var(months[start_ind:end_ind],
                     rmses[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("month")
        plt.ylabel("{} RMSE ({})".format(name, units[name]))
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        month_weights = [m / 365. for m in month_days]
        print(name, " has reference mean: ", sum([ref_means[name][j] * month_weights[j] for j in range(START_AVG_MONTH-START_MONTH, END_AVG_MONTH-START_MONTH+1)]),
              file=log_file)
        for i in range(case_num):
            print(name, " has case ", TEST_CASES[i].short_name, " mean: ", sum([test_means[name][i,j] * month_weights[j] for j in range(START_AVG_MONTH-START_MONTH, END_AVG_MONTH-START_MONTH+1)]),
                  file=log_file)
            print(name, " has difference mean: ", sum([diff_means[name][i,j] * month_weights[j] for j in range(START_AVG_MONTH-START_MONTH, END_AVG_MONTH-START_MONTH+1)]),
                  file=log_file)

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
}
names = list(units.keys())
scales = dict()
for name in names:
    scales[name] = 1.
scales['SWCF'] = -1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.
scales['PRECT'] = 1000.*86400.

log_plot_names = []#'AODABS', 'AODVIS', 'AODUV']

plot_vars_over_time(names, units, scales, log_plot_names)

log_file.close()
