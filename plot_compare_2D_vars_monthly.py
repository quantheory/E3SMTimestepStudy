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

START_YEAR = 1
START_MONTH = 3
END_YEAR = 4
END_MONTH = 2

MONTHLY_FILE_LOC="/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon"

USE_PRESAER=False

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
if USE_PRESAER:
    suffix += '_presaer'

log_file = open("plot_2D_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", MONTHLY_FILE_LOC, START_MONTH, END_MONTH)
    TEST_CASES = [
    ]
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", MONTHLY_FILE_LOC, START_MONTH, END_MONTH)
    TEST_CASES = [
        E3SMCaseOutput("timestep_all_10s", "ALL10", MONTHLY_FILE_LOC, START_MONTH, END_MONTH),
    ]

rfile0 = nc4.Dataset(REF_CASE.get_monthly_file_name(START_MONTH, START_YEAR), 'r')
lat = rfile0['lat'][:]
lon = rfile0['lon'][:]
nlat = len(lat)
nlon = len(lon)
rfile0.close()

month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_weights = [month_days[months[imonth] - 1] for imonth in imonths]
weight_norm = sum(month_weights)
month_weights = [weight / weight_norm for weight in month_weights]

def get_overall_averages(ref_case, test_cases, months, varnames, scales):
    case_num = len(test_cases)
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
            ref_means[name] += month_weights[imonth] * ref_monthly[name]
            for icase in range(case_num):
                test_means[icase][name] += month_weights[imonth] * test_monthly[icase][name]
                diff_means[icase][name] += month_weights[imonth] * diff_monthly[icase][name]

    for name in varnames:
        ref_means[name] *= scales[name]
        for icase in range(case_num):
            test_means[icase][name] *= scales[name]
            diff_means[icase][name] *= scales[name]

    return (ref_means, test_means, diff_means)

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
    'U10': r'$m/s$'
}
varnames = list(units.keys())
scales = dict()
for name in varnames:
    scales[name] = 1.
scales['SWCF'] = -1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.
scales['PRECT'] = 1000.*86400.

diff_lims = {
    'OMEGA500': 0.05,
    'TAU': 0.1,
    'PRECT': 4.,
}

ref_means, test_means, diff_means = get_overall_averages(REF_CASE, TEST_CASES, months, varnames, scales)

for name in varnames:
    plot_name = name
    if name in plot_names:
        plot_name = plot_names[name]
    clim_val = [ref_means[name].min(), ref_means[name].max()]
    clim_diff = 0.
    for icase in range(len(TEST_CASES)):
        clim_val[0] = min(clim_val[0], test_means[icase][name].min())
        clim_val[1] = max(clim_val[1], test_means[icase][name].max())
        clim_diff = max(clim_diff, -diff_means[icase][name].min())
        clim_diff = max(clim_diff, diff_means[icase][name].max())

    if name in diff_lims:
        clim_diff = diff_lims[name]

    plt.pcolormesh(lon[:], lat[:], ref_means[name])
    bmap.drawcoastlines()
    ax = plt.gca()
    ax.set_xticks([0., 90., 180., 270., 360.])
    ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
    ax.set_yticks([60., 30., 0., -30., -60.])
    ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
    plt.axis('tight')
    plt.xlim([0., 360.])
    plt.colorbar()
    plt.clim(clim_val[0], clim_val[1])
    plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, REF_CASE.short_name, units[name],
                                                                  day_str(START_MONTH), day_str(START_YEAR),
                                                                  day_str(END_MONTH), day_str(END_YEAR)))
    plt.savefig('{}_{}{}.png'.format(name, REF_CASE.short_name, suffix))
    plt.close()

    for icase in range(len(TEST_CASES)):
        case_name = TEST_CASES[icase].short_name
        plt.pcolormesh(lon[:], lat[:], test_means[icase][name])
        bmap.drawcoastlines()
        ax = plt.gca()
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yticks([60., 30., 0., -30., -60.])
        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.axis('tight')
        plt.xlim([0., 360.])
        plt.colorbar()
        plt.clim(clim_val[0], clim_val[1])
        plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, case_name, units[name],
                                                                      day_str(START_MONTH), day_str(START_YEAR),
                                                                      day_str(END_MONTH), day_str(END_YEAR)))
        plt.savefig('{}_{}{}.png'.format(name, case_name, suffix))
        plt.close()

        plt.pcolormesh(lon[:], lat[:], diff_means[icase][name], cmap=cmap)
        bmap.drawcoastlines()
        ax = plt.gca()
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yticks([60., 30., 0., -30., -60.])
        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.axis('tight')
        plt.xlim([0., 360.])
        plt.colorbar()
        plt.clim(-clim_diff, clim_diff)
        plt.title("Mean difference in {}\nfor case {} ({}, months {}/{} - {}/{})".format(plot_name, case_name, units[name],
                                                                                         day_str(START_MONTH), day_str(START_YEAR),
                                                                                         day_str(END_MONTH), day_str(END_YEAR)))
        plt.savefig('{}_diff_{}{}.png'.format(name, case_name, suffix))
        plt.close()
