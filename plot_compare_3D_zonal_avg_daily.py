#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

cmap = plt.get_cmap('coolwarm')

def forward(a):
    a = np.deg2rad(a)
    return np.sin(a)

def inverse(a):
    a = np.arcsin(a)
    return np.rad2deg(a)

START_DAY = 3
END_DAY = 15

DAILY_FILE_LOC="/p/lscratchh/santos36/timestep_daily_avgs_lat_lon"

USE_PRESAER = False

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)

suffix = '_d{}-{}'.format(day_str(START_DAY), day_str(END_DAY))

suffix += '_zonal'

log_file = open("plot_zonal_log{}.txt".format(suffix), 'w')

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
        E3SMCaseOutput("timestep_CLUBB_MG2_Strang", "CLUBBMICROSTR", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_Strang_60s", "CLUBBMICROSTR60", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_ZM_10s", "ZM10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_ZM_300s", "ZM300", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_rad_10s", "ALLRAD10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_300s", "ALL300", DAILY_FILE_LOC, START_DAY, END_DAY),
    ]

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
lat = rfile0['lat'][:]
lon = rfile0['lon'][:]
ilev = rfile0['ilev'][:]
nlat = len(lat)
nlon = len(lon)
nlev = len(ilev) - 1
rfile0.close()

def get_overall_averages(ref_case, test_cases, days, varnames, scales):
    case_num = len(test_cases)
    day_num = len(days)
    ref_means = dict()
    for name in varnames:
        ref_means[name] = np.zeros((nlev, nlat, nlon))
    test_means = []
    diff_means = []
    for case in test_cases:
        next_test_means = dict()
        next_diff_means = dict()
        for name in varnames:
            next_test_means[name] = np.zeros((nlev, nlat, nlon))
            next_diff_means[name] = np.zeros((nlev, nlat, nlon))
        test_means.append(next_test_means)
        diff_means.append(next_diff_means)

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

    for day in days:
        ref_daily, test_daily, diff_daily = ref_case.compare_daily_averages(test_cases, day, varnames_read)
        if "PRECT" in varnames:
            ref_daily["PRECT"] = ref_daily["PRECL"] + ref_daily["PRECC"]
            for icase in range(case_num):
                test_daily[icase]["PRECT"] = test_daily[icase]["PRECL"] + test_daily[icase]["PRECC"]
                diff_daily[icase]["PRECT"] = diff_daily[icase]["PRECL"] + diff_daily[icase]["PRECC"]
        if "TAU" in varnames:
            ref_daily["TAU"] = np.sqrt(ref_daily["TAUX"]**2 + ref_daily["TAUY"]**2)
            for icase in range(case_num):
                test_daily[icase]["TAU"] = np.sqrt(test_daily[icase]["TAUX"]**2 + test_daily[icase]["TAUY"]**2)
                diff_daily[icase]["TAU"] = test_daily[icase]["TAU"] - ref_daily["TAU"]
        for name in varnames:
            for jlev in range(nlev):
                ref_means[name][jlev,:,:] += ref_daily[name][jlev,:,:]
                for icase in range(case_num):
                    test_means[icase][name][jlev,:,:] += test_daily[icase][name][jlev,:,:]
                    diff_means[icase][name][jlev,:,:] += diff_daily[icase][name][jlev,:,:]

    for name in varnames:
        ref_means[name] *= scales[name]/day_num
        for icase in range(case_num):
            test_means[icase][name] *= scales[name]/day_num
            diff_means[icase][name] *= scales[name]/day_num

    return (ref_means, test_means, diff_means)


plot_names = {
    'AQRAIN': 'rain mixing ratio',
    'AQSNOW': 'snow mixing ratio',
    'AREI': 'ice effective radius',
    'AREL': 'droplet effective radius',
    'CLDICE': 'cloud ice mixing ratio',
    'CLDLIQ': 'cloud liquid mixing ratio',
    'CLOUD': "cloud fraction",
    'Q': "specific humidity",
    'QRL': 'longwave heating rate',
    'QRS': 'shortwave heating rate',
    'RELHUM': "relative humidity",
    'T': "temperature",
    'U': "zonal wind",
    'V': "meridional wind",
}

units = {
    'AQRAIN': r'$mg/kg$',
    'AQSNOW': r'$mg/kg$',
    'AREI': r'micron',
    'AREL': r'micron',
    'CLDICE': r'$g/kg$',
    'CLDLIQ': r'$g/kg$',
    'CLOUD': r'fraction',
    'Q': r'$g/kg$',
    'QRL': r'$K/day$',
    'QRS': r'$K/day$',
    'RELHUM': r'%',
    'T': r'$K$',
    'U': r'$m/s$',
    'V': r'$m/s$',
}
varnames = list(units.keys())
scales = dict()
for name in varnames:
    scales[name] = 1.
scales['AQRAIN'] = 1.e6
scales['AQSNOW'] = 1.e6
scales['CLDICE'] = 1000.
scales['CLDLIQ'] = 1000.
scales['Q'] = 1000.
scales['QRL'] = 86400.
scales['QRS'] = 86400.

diff_lims = {
    'AQSNOW': 4.,
    'AREI': 5.,
    'AREL': 2.,
    'CLDLIQ': 0.01,
    'CLOUD': 0.1,
    'QRL': 2.,
    'QRS': 0.3,
    'RELHUM': 15.,
    'T': 2.,
}

PLOT_TOP = 100.
itop = 0
for level in ilev:
    if level > PLOT_TOP:
        break
    itop += 1

def zonal_average(x):
    return np.mean(x[itop:,:,:], axis=2)

ref_means, test_means, diff_means = get_overall_averages(REF_CASE, TEST_CASES, days, varnames, scales)

plot_ilev = ilev[itop:]

for name in varnames:
    plot_name = name
    if name in plot_names:
        plot_name = plot_names[name]

    ref_plot_var = zonal_average(ref_means[name])
    clim_val = [ref_plot_var.min(), ref_plot_var.max()]
    clim_diff = 0.
    for icase in range(len(TEST_CASES)):
        test_plot_var = zonal_average(test_means[icase][name])
        diff_plot_var = zonal_average(diff_means[icase][name])
        clim_val[0] = min(clim_val[0], test_plot_var.min())
        clim_val[1] = max(clim_val[1], test_plot_var.max())
        clim_diff = max(clim_diff, - diff_plot_var.min())
        clim_diff = max(clim_diff, diff_plot_var.max())

    if name in diff_lims:
        clim_diff = diff_lims[name]

    plt.pcolor(lat[1:], plot_ilev, ref_plot_var[:,1:-1])
    ax = plt.gca()
    ylim = ax.get_ylim()
    ax.set_ylim([ylim[1], ylim[0]])
    plt.ylabel("Pressure (mb)")
    ax.set_xscale('function', functions=(forward, inverse))
    ax.set_xticks([60., 30., 15., 0., -15., -30., -60.])
    ax.set_xticklabels(['60N', '30N', '15N', '0', '15S', '30S', '60S'])
    plt.axis('tight')
    plt.colorbar()
    plt.clim(clim_val[0], clim_val[1])
    plt.title("{} for case {}\n({}, days {}-{})".format(plot_name, REF_CASE.short_name, units[name],
                                                        START_DAY, END_DAY))
    plt.savefig('{}_{}{}.png'.format(name, REF_CASE.short_name, suffix))
    plt.close()

    for icase in range(len(TEST_CASES)):
        test_plot_var = zonal_average(test_means[icase][name])
        diff_plot_var = zonal_average(diff_means[icase][name])
        case_name = TEST_CASES[icase].short_name

        plt.pcolor(lat[1:], plot_ilev, test_plot_var[:,1:-1])
        ax = plt.gca()
        ylim = ax.get_ylim()
        ax.set_ylim([ylim[1], ylim[0]])
        plt.ylabel("Pressure (mb)")
        ax.set_xscale('function', functions=(forward, inverse))
        ax.set_xticks([60., 30., 15., 0., -15., -30., -60.])
        ax.set_xticklabels(['60N', '30N', '15N', '0', '15S', '30S', '60S'])
        plt.axis('tight')
        plt.colorbar()
        plt.clim(clim_val[0], clim_val[1])
        plt.title("{} for case {}\n({}, days {}-{})".format(plot_name, case_name, units[name],
                                                                      START_DAY, END_DAY))
        plt.savefig('{}_{}{}.png'.format(name, case_name, suffix))
        plt.close()

        plt.pcolor(lat[1:], plot_ilev, diff_plot_var[:,1:-1], cmap=cmap)
        ax = plt.gca()
        ylim = ax.get_ylim()
        ax.set_ylim([ylim[1], ylim[0]])
        plt.ylabel("Pressure (mb)")
        ax.set_xscale('function', functions=(forward, inverse))
        ax.set_xticks([60., 30., 15., 0., -15., -30., -60.])
        ax.set_xticklabels(['60N', '30N', '15N', '0', '15S', '30S', '60S'])
        plt.axis('tight')
        plt.colorbar()
        plt.clim(-clim_diff, clim_diff)
        plt.title("Mean difference in {}\nfor case {} ({}, days {}-{})".format(plot_name, case_name, units[name],
                                                                               START_DAY, END_DAY))
        plt.savefig('{}_diff_{}{}.png'.format(name, case_name, suffix))
        plt.close()
