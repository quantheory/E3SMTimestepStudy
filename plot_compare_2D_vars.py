#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

cmap = plt.get_cmap('coolwarm')
bmap = basemap.Basemap(lon_0=180.)

START_DAY = 3
END_DAY = 15

DAILY_FILE_LOC="/p/lscratchh/santos36/timestep_daily_avgs_lat_lon"

USE_PRESAER=True

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)

suffix = '_d{}-{}'.format(day_str(START_DAY), day_str(END_DAY))
if USE_PRESAER:
    suffix += '_presaer'

log_file = open("plot_2D_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", DAILY_FILE_LOC)
    TEST_CASES = [
        E3SMCaseOutput("timestep_presaer_all_10s", "ALL10PA", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s", "CLUBBMICRO10PA", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_presaer_ZM_10s", "ZM10PA", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s_ZM_10s", "CLUBBMICRO10ZM10PA", DAILY_FILE_LOC),
    ]
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC)
    TEST_CASES = [
        E3SMCaseOutput("timestep_all_10s", "ALL10", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_dyn_10s", "DYN10", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_MG2_10s", "MICRO10", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_CLUBB_10s", "CLUBB10", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_CLUBB_10s_MG2_10s", "CLUBB10MICRO10", DAILY_FILE_LOC),
        E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_ZM_10s", "ZM10", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_CLUBB_MG2_10s_ZM_10s", "CLUBBMICRO10ZM10", DAILY_FILE_LOC),
        #E3SMCaseOutput("timestep_ZM_300s", "ZM300", DAILY_FILE_LOC),
    ]

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
lat = rfile0['lat'][:]
lon = rfile0['lon'][:]
nlat = len(lat)
nlon = len(lon)
rfile0.close()

def get_overall_averages(ref_case, test_cases, days, varnames):
    case_num = len(test_cases)
    day_num = len(days)
    ref_means = dict()
    for name in varnames:
        ref_means[name] = np.zeros((nlat, nlon))
    test_means = []
    diff_means = []
    for case in test_cases:
        next_test_means = dict()
        next_diff_means = dict()
        for name in varnames:
            next_test_means[name] = np.zeros((nlat, nlon))
            next_diff_means[name] = np.zeros((nlat, nlon))
        test_means.append(next_test_means)
        diff_means.append(next_diff_means)

    for day in days:
        ref_daily, test_daily, diff_daily = ref_case.compare_daily_averages(test_cases, day, varnames)
        for name in varnames:
            ref_means[name] += ref_daily[name]
            for icase in range(case_num):
                test_means[icase][name] += test_daily[icase][name]
                diff_means[icase][name] += diff_daily[icase][name]

    for name in varnames:
        ref_means[name] /= day_num
        for icase in range(case_num):
            test_means[icase][name] /= day_num
            diff_means[icase][name] /= day_num

    return (ref_means, test_means, diff_means)

units = {
    'LWCF': r'$W/m^2$',
    'SWCF': r'$W/m^2$',
    'PRECC': r'$mm/day$',
    'PRECL': r'$mm/day$',
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
}
varnames = list(units.keys())
scales = dict()
for name in varnames:
    scales[name] = 1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.

ref_means, test_means, diff_means = get_overall_averages(REF_CASE, TEST_CASES, days, varnames)

for name in varnames:
    clim_val = [ref_means[name].min(), ref_means[name].max()]
    clim_diff = [1.e36, 1.e-36] # Kind of hacky...
    for icase in range(len(TEST_CASES)):
        clim_val[0] = min(clim_val[0], test_means[icase][name].min())
        clim_val[1] = max(clim_val[1], test_means[icase][name].max())
        clim_diff[0] = min(clim_diff[0], diff_means[icase][name].min())
        clim_diff[1] = max(clim_diff[1], diff_means[icase][name].max())

    plt.pcolormesh(lon[:], lat[:], ref_means[name])
    bmap.drawcoastlines()
    ax = plt.gca()
    plt.axis('tight')
    ax.set_xticks([0., 90., 180., 270., 360.])
    ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
    ax.set_yticks([60., 30., 0., -30., -60.])
    ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
    plt.colorbar()
    plt.clim(clim_val[0], clim_val[1])
    plt.title("{} for case {} (days {}-{})".format(name, REF_CASE.short_name, START_DAY, END_DAY))
    plt.savefig('{}_{}{}.png'.format(name, REF_CASE.short_name, suffix))
    plt.close()

    for icase in range(len(TEST_CASES)):
        case_name = TEST_CASES[icase].short_name
        plt.pcolormesh(lon[:], lat[:], test_means[icase][name])
        bmap.drawcoastlines()
        ax = plt.gca()
        plt.axis('tight')
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yticks([60., 30., 0., -30., -60.])
        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.colorbar()
        plt.clim(clim_val[0], clim_val[1])
        plt.title("{} for case {} (days {}-{})".format(name, case_name, START_DAY, END_DAY))
        plt.savefig('{}_{}{}.png'.format(name, case_name, suffix))
        plt.close()

        plt.pcolormesh(lon[:], lat[:], diff_means[icase][name], cmap=cmap)
        bmap.drawcoastlines()
        ax = plt.gca()
        plt.axis('tight')
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yticks([60., 30., 0., -30., -60.])
        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.colorbar()
        plt.clim(clim_diff[0], clim_diff[1])
        plt.title("Mean difference in {} for case {} (days {}-{})".format(name, case_name, START_DAY, END_DAY))
        plt.savefig('{}_diff_{}{}.png'.format(name, case_name, suffix))
        plt.close()
