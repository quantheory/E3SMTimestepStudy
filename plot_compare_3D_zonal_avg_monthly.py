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

START_YEAR = 1
START_MONTH = 3
END_YEAR = 4
END_MONTH = 2

MONTHLY_FILE_LOC="/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon"

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


suffix = '_y{}m{}-y{}m{}'.format(day_str(START_YEAR),
                                 day_str(START_MONTH),
                                 day_str(END_YEAR),
                                 day_str(END_MONTH))

suffix += '_zonal'

log_file = open("plot_zonal_log{}.txt".format(suffix), 'w')

REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", MONTHLY_FILE_LOC, START_MONTH, END_MONTH)
TEST_CASES = [
    E3SMCaseOutput("timestep_all_10s", "ALL10", MONTHLY_FILE_LOC, START_MONTH, END_MONTH),
]

rfile0 = nc4.Dataset(REF_CASE.get_monthly_file_name(START_MONTH, START_YEAR), 'r')
lat = rfile0['lat'][:]
lon = rfile0['lon'][:]
ilev = rfile0['ilev'][:]
nlat = len(lat)
nlon = len(lon)
nlev = len(ilev) - 1
rfile0.close()

month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_weights = [month_days[months[imonth] - 1] for imonth in imonths]
weight_norm = sum(month_weights)
month_weights = [weight / weight_norm for weight in month_weights]

def get_overall_averages(ref_case, test_cases, months, varnames, scales):
    case_num = len(test_cases)
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

    for imonth in imonths:
        month = months[imonth]
        year = years[imonth]
        ref_monthly, test_monthly, diff_monthly = ref_case.compare_monthly_averages(test_cases, month, year, varnames_read)
        if "PRECT" in varnames:
            ref_monthly["PRECT"] = ref_monthly["PRECL"] + ref_monthly["PRECC"]
            for icase in range(case_num):
                test_monthly[icase]["PRECT"] = test_monthly[icase]["PRECL"] + test_monthly[icase]["PRECC"]
                diff_monthly[icase]["PRECT"] = diff_monthly[icase]["PRECL"] + diff_monthly[icase]["PRECC"]
        if "TAU" in varnames:
            ref_monthly["TAU"] = np.sqrt(ref_monthly["TAUX"]**2 + ref_monthly["TAUY"]**2)
            for icase in range(case_num):
                test_monthly[icase]["TAU"] = np.sqrt(test_monthly[icase]["TAUX"]**2 + test_monthly[icase]["TAUY"]**2)
                diff_monthly[icase]["TAU"] = test_monthly[icase]["TAU"] - ref_monthly["TAU"]
        for name in varnames:
            for jlev in range(nlev):
                ref_means[name][jlev,:,:] += month_weights[imonth] * ref_monthly[name][jlev,:,:]
                for icase in range(case_num):
                    test_means[icase][name][jlev,:,:] += month_weights[imonth] * test_monthly[icase][name][jlev,:,:]
                    diff_means[icase][name][jlev,:,:] += month_weights[imonth] * diff_monthly[icase][name][jlev,:,:]

    for name in varnames:
        ref_means[name] *= scales[name]
        for icase in range(case_num):
            test_means[icase][name] *= scales[name]
            diff_means[icase][name] *= scales[name]

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
    'U': r'$m/s$'
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

PLOT_TOP = 100.
itop = 0
for level in ilev:
    if level > PLOT_TOP:
        break
    itop += 1

def zonal_average(x):
    return np.mean(x[itop:,:,:], axis=2)

ref_means, test_means, diff_means = get_overall_averages(REF_CASE, TEST_CASES, months, varnames, scales)

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
    plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, REF_CASE.short_name, units[name],
                                                                  day_str(START_MONTH), day_str(START_YEAR),
                                                                  day_str(END_MONTH), day_str(END_YEAR)))
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
        plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, case_name, units[name],
                                                                      day_str(START_MONTH), day_str(START_YEAR),
                                                                      day_str(END_MONTH), day_str(END_YEAR)))
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
        plt.title("Mean difference in {}\nfor case {} ({}, months {}/{} - {}/{})".format(plot_name, case_name, units[name],
                                                                                         day_str(START_MONTH), day_str(START_YEAR),
                                                                                         day_str(END_MONTH), day_str(END_YEAR)))
        plt.savefig('{}_diff_{}{}.png'.format(name, case_name, suffix))
        plt.close()
