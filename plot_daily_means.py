#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

start_day = 1
end_day = 30

days = list(range(start_day, end_day+1))
ndays = len(days)

dayform = "{:02d}"
day_strs = [dayform.format(day) for day in days]
day_strs.append(dayform.format(end_day+1))

suffix = '_d{}-{}'.format(day_strs[0], day_strs[-2])

log_file = open("plot_daily_log{}.txt".format(suffix), 'w')

times = ["03600", "07200", "10800", "14400", "18000", "21600", "25200", "28800",
         "32400", "36000", "39600", "43200", "46800", "50400", "54000", "57600",
         "61200", "64800", "68400", "72000", "75600", "79200", "82800", "00000"]

ntimes = len(times)

REF_CASE_NAME="timestep_ctrl"
REF_CASE_LABEL="CTRL"
REF_CASE_RUNDIR="/p/lscratchh/santos36/ACME/{}/run/".format(REF_CASE_NAME)
TEST_CASE_LABELS={
    "timestep_all_10s": "ALL10",
    "timestep_dyn_10s": "DYN10",
    "timestep_MG2_10s": "MICRO10",
    "timestep_CLUBB_10s": "CLUBB10",
    "timestep_CLUBB_MG2_10s": "CLUBBMICRO10",
    "timestep_CLUBB_MG2_60s": "CLUBBMICRO60",
#    "timestep_ZM_10s": "ZM10",
    "timestep_ZM_300s": "ZM300",
}
TEST_CASE_NAMES=list(TEST_CASE_LABELS.keys())
TEST_CASE_RUNDIRS=["/p/lscratchh/santos36/ACME/{}/run/".format(name)
                   for name in TEST_CASE_NAMES]

case_num = len(TEST_CASE_NAMES)

REF_FILE_NAMES = []
TEST_FILE_NAMES = [[] for i in range(case_num)]

for iday in range(ndays):
    REF_FILE_NAMES.append([])
    for i in range(case_num):
        TEST_FILE_NAMES[i].append([])
    for time in times:
        if time != "00000":
            day_str = day_strs[iday]
        else:
            day_str = day_strs[iday+1]
        REF_FILE_NAMES[iday].append('{}/{}.cam.h0.0001-01-{}-{}.nc'.format(REF_CASE_RUNDIR, REF_CASE_NAME, day_str, time))
        for i in range(case_num):
            TEST_FILE_NAMES[i][iday].append('{}/{}.cam.h0.0001-01-{}-{}.nc'.format(TEST_CASE_RUNDIRS[i], TEST_CASE_NAMES[i], day_str, time))

rfile0 = nc4.Dataset(REF_FILE_NAMES[0][0], 'r')
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
area_sum = area.sum()
weights = area/area_sum
rfile0.close()

def average_2D_var_over_files(file_names, varnames):
    var_avg = dict()
    for varname in varnames:
        var_avg[varname] = np.zeros((ncol,))
    for file_name in file_names:
        ncfile = nc4.Dataset(file_name, 'r')
        for varname in varnames:
            if "missing_value" in ncfile[varname].ncattrs():
                file_var = ncfile[varname][0,:]
                miss = ncfile[varname].missing_value
                # Following appears not to work due to rounding errors.
                # var_avg[varname] += np.where(file_var == miss, 0., file_var)
                # Kludge for optical depths.
                var_avg[varname] += np.where(file_var > 1.e35, 0., file_var)
            else:
                var_avg[varname] += ncfile[varname][0,:]
        ncfile.close()
    for varname in varnames:
        var_avg[varname] /= len(file_names)
    return var_avg

def calc_2D_var_stats(ref_names, test_names, varnames):
    ref_time_avg = average_2D_var_over_files(ref_names, varnames)
    test_time_avgs = []
    diff_time_avgs = []
    for i in range(case_num):
        test_time_avgs.append(average_2D_var_over_files(test_names[i], varnames))
        next_diff_time = dict()
        for varname in varnames:
            next_diff_time[varname] = test_time_avgs[i][varname] - ref_time_avg[varname]
        diff_time_avgs.append(next_diff_time)
    ref_avg = dict()
    test_avgs = dict()
    diff_avgs = dict()
    rmses = dict()
    for varname in varnames:
        ref_avg[varname] = (ref_time_avg[varname] * weights).sum()
        test_avgs[varname] = []
        diff_avgs[varname] = []
        rmses[varname] = []
        for i in range(case_num):
            test_avgs[varname].append((test_time_avgs[i][varname] * weights).sum())
            diff_avgs[varname].append((diff_time_avgs[i][varname] * weights).sum())
            #assert np.isclose(diff_avgs[varname][i], test_avgs[varname][i] - ref_avg[varname]), \
            #    "Problem with diff of variable {} from case {}".format(varname, TEST_CASE_NAMES[i])
            rmses[varname].append(np.sqrt((diff_time_avgs[i][varname]**2 * weights).sum()))
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
        print("On day: ", iday+start_day, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_2D_var_stats(REF_FILE_NAMES[iday],
                                                                                   [cfiles[iday] for cfiles in TEST_FILE_NAMES],
                                                                                   names)
        for name in names:
            ref_means[name][iday] = ref_mean[name]*scales[name]
            for i in range(case_num):
                test_means[name][i,iday] = test_case_means[name][i]*scales[name]
                diff_means[name][i,iday] = diff_case_means[name][i]*scales[name]
                rmses[name][i,iday] = case_rmses[name][i]*scales[name]

    for name in names:
        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        plot_var(days, ref_means[name], label=REF_CASE_LABEL)
        for i in range(case_num):
            plot_var(days, test_means[name][i,:],
                     label=TEST_CASE_LABELS[TEST_CASE_NAMES[i]])
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} ({})".format(name, units[name]))
        plt.legend()
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            plot_var(days, diff_means[name][i,:],
                     label=TEST_CASE_LABELS[TEST_CASE_NAMES[i]])
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} difference ({})".format(name, units[name]))
        plt.legend()
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            plt.plot(days, rmses[name][i,:],
                     label=TEST_CASE_LABELS[TEST_CASE_NAMES[i]])
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("{} RMSE ({})".format(name, units))
        plt.legend()
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        print(name, " has reference mean: ", sum(ref_means[name])/ndays,
              file=log_file)
        for i in range(case_num):
            print(name, " has case ", TEST_CASE_NAMES[i], " mean: ", sum(test_means[name][i])/ndays,
                  file=log_file)
            print(name, " has difference mean: ", sum(diff_means[name][i])/ndays,
                  file=log_file)

units = {
    'LWCF': r'$W/m^2$',
    'SWCF': r'$W/m^2$',
    'PRECC': r'$mm/day$',
    'PRECL': r'$mm/day$',
    'TGCLDCWP': r'$kg/m^2$',
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
}
names = list(units.keys())
scales = dict()
for name in names:
    scales[name] = 1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.

log_plot_names = []#'AODABS', 'AODVIS', 'AODUV']

plot_vars_over_time(names, units, scales, log_plot_names)

log_file.close()
