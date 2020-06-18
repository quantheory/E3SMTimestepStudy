#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

START_DAY = 1
END_DAY = 15
END_ZM10S_DAY = 19

START_AVG_DAY = 3
END_AVG_DAY = 15

DAILY_FILE_LOC="/p/lscratchh/santos36/timestep_daily_avgs/"

USE_PRESAER=False
TROPICS_ONLY=False

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)
navgdays = END_AVG_DAY - START_AVG_DAY + 1

suffix = '_d{}-{}'.format(day_str(START_DAY), day_str(END_DAY))
if USE_PRESAER:
    suffix += '_presaer'
if TROPICS_ONLY:
    suffix += '_tropics'

log_file = open("plot_daily_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_presaer_all_10s", "ALL10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s", "CLUBBMICRO10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_ZM_10s", "ZM10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s_ZM_10s", "CLUBBMICRO10ZM10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_cld_10s", "CLD10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
    ]
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_all_10s", "ALL10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_dyn_10s", "DYN10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_MG2_10s", "MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_10s", "CLUBB10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_10s_MG2_10s", "CLUBB10MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_ZM_10s", "ZM10", DAILY_FILE_LOC, START_DAY, END_ZM10S_DAY),
        E3SMCaseOutput("timestep_ZM_300s", "ZM300", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_rad_10s", "ALLRAD10", DAILY_FILE_LOC, START_DAY, END_DAY),
    ]

case_num = len(TEST_CASES)

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
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

def calc_2D_var_stats(ref_case, test_cases, day, varnames):
    varnames_read = [name for name in varnames if name != "PRECT"]
    if "PRECT" in varnames:
        if "PRECL" not in varnames:
            varnames_read.append("PRECL")
        if "PRECC" not in varnames:
            varnames_read.append("PRECC")
    ref_time_avg, test_time_avgs, diff_time_avgs = ref_case.compare_daily_averages(test_cases, day, varnames_read)
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
            if test_cases[i].day_is_available(day):
                test_avgs[varname].append((test_time_avgs[i][varname] * weights).sum())
                diff_avgs[varname].append((diff_time_avgs[i][varname] * weights).sum())
                assert np.isclose(diff_avgs[varname][i], test_avgs[varname][i] - ref_avg[varname]), \
                    "Problem with diff of variable {} from case {}".format(varname, TEST_CASE_NAMES[i])
                rmses[varname].append(np.sqrt((diff_time_avgs[i][varname]**2 * weights).sum()))
            else:
                test_avgs[varname].append(None)
                diff_avgs[varname].append(None)
                rmses[varname].append(None)
    return (ref_avg, test_avgs, diff_avgs, rmses)

def plot_vars_over_time(names, units, scales, log_plot_names):
    ref_means = dict()
    test_means = dict()
    diff_means = dict()
    rmses = dict()
    for name in names:
        ref_means[name] = np.zeros((ndays,))
        test_means[name] = np.zeros((case_num, ndays))
        diff_means[name] = np.zeros((case_num, ndays))
        rmses[name] = np.zeros((case_num, ndays))

    for iday in range(ndays):
        day = days[iday]
        print("On day: ", day, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_2D_var_stats(REF_CASE, TEST_CASES, day, names)
        for name in names:
            ref_means[name][iday] = ref_mean[name]*scales[name]
            for i in range(case_num):
                if TEST_CASES[i].day_is_available(day):
                    test_means[name][i,iday] = test_case_means[name][i]*scales[name]
                    diff_means[name][i,iday] = diff_case_means[name][i]*scales[name]
                    rmses[name][i,iday] = case_rmses[name][i]*scales[name]

    for name in names:
        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        plot_var(days, ref_means[name], label=REF_CASE.short_name)
        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     test_means[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} ({})".format(name, units[name]))
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     diff_means[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} difference ({})".format(name, units[name]))
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     rmses[name][i,start_ind:end_ind],
                     label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("{} RMSE ({})".format(name, units[name]))
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        print(name, " has reference mean: ", sum(ref_means[name][START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
              file=log_file)
        for i in range(case_num):
            print(name, " has case ", TEST_CASES[i].short_name, " mean: ", sum(test_means[name][i,START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
                  file=log_file)
            print(name, " has difference mean: ", sum(diff_means[name][i,START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
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
