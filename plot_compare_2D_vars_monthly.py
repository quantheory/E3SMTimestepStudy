#!/usr/bin/env python

from os.path import join

from functools import partial

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

cmap = plt.get_cmap('coolwarm')
bmap = basemap.Basemap(lon_0=180.)

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

MONTHLY_FILE_LOC = "/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon"

USE_PRESAER = False

PLOT_FILE_TYPE = "png"

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
lev = rfile0['lev'][:]
nlat = len(lat)
nlon = len(lon)
nlev = len(lev)
rfile0.close()

month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_weights = [month_days[months[imonth] - 1] for imonth in imonths]
weight_norm = sum(month_weights)
month_weights = [weight / weight_norm for weight in month_weights]

def get_overall_averages(ref_case, test_cases, months, varnames, scales):
    case_num = len(test_cases)
    ref_means = dict()
    for name in varnames:
        if name in vars_3D:
            ref_means[name] = np.zeros((nlev, nlat, nlon))
        else:
            ref_means[name] = np.zeros((nlat, nlon))
    test_means = []
    diff_means = []
    for case in test_cases:
        next_test_means = dict()
        next_diff_means = dict()
        for name in varnames:
            if name in vars_3D:
                next_test_means[name] = np.zeros((nlev, nlat, nlon))
                next_diff_means[name] = np.zeros((nlev, nlat, nlon))
            else:
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
            if name in vars_3D:
                for ilev in range(nlev):
                    ref_means[name][ilev,:,:] += month_weights[imonth] * ref_monthly[name][ilev,:,:]
                    for icase in range(case_num):
                        test_means[icase][name][ilev,:,:] += month_weights[imonth] * test_monthly[icase][name][ilev,:,:]
                        diff_means[icase][name][ilev,:,:] += month_weights[imonth] * diff_monthly[icase][name][ilev,:,:]
            else:
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
    'PRECL': "large-scale precipitation",
    'PRECE': "extreme precipitation",
    'PRECT': "total precipitation",
    'TGCLDIWP': "ice water path",
    'TGCLDLWP': "liquid water path",
    'CLDTOT': "cloud area fraction",
    'CLDLOW': "low cloud area fraction",
    'CLDMED': "mid-level cloud area fraction",
    'CLDHGH': "high cloud area fraction",
    'LHFLX': "latent heat flux",
    'SHFLX': "sensible heat flux",
    'TAU': "surface wind stress",
    'TS': "surface temperature",
    'PSL': "sea level pressure",
    'OMEGA500': "vertical velocity at 500 mb",
    'U10': "10 meter wind speed",
    'RELHUM': "surface relative humidity",
    'Q': "specific humidity",
    'CLDLIQ': "lowest level cloud liquid",
    'T': "lowest level temperature",
    'CLOUD': "lowest level cloud fraction",
    'TMQ': "precipitable water",
}

units = {
    'LWCF': r'$W/m^2$',
    'SWCF': r'$W/m^2$',
    'PRECC': r'$mm/day$',
    'PRECL': r'$mm/day$',
    'PRECE': r'$mm/day$',
    'PRECT': r'$mm/day$',
    'TGCLDIWP': r'$g/m^2$',
    'TGCLDLWP': r'$g/m^2$',
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
    'CLDLIQ': r"$g/kg$",
    'T': r'$K$',
    'CLOUD': r'$fraction$',
    'TMQ': r'$kg/m^2$',
}
varnames = list(units.keys())
# Only can plot extreme precipitation for a subset of months.
varnames.remove('PRECE')
scales = dict()
for name in varnames:
    scales[name] = 1.
scales['TGCLDIWP'] = 1000.
scales['TGCLDLWP'] = 1000.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.
scales['PRECE'] = 1000.*86400.
scales['PRECT'] = 1000.*86400.
scales['Q'] = 1000.
scales['CLDLIQ'] = 1000.

diff_lims = {
    'OMEGA500': 0.05,
    'TAU': 0.1,
    'PRECC': 4.,
    'PRECL': 4.,
    'PRECT': 4.,
    'PSL': 200.
}

vars_3D = [
    'RELHUM',
    'Q',
    'CLDLIQ',
    'T',
    'CLOUD',
]

PRECIP_OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip_lat_lon/"

out_file_template = "{}.freq.00{}-{}.nc"

# Threshold for precipitation to be considered "extreme", in mm/day.
PRECE_THRESHOLD = 97.

def get_prece(case_name):
    file_name = out_file_template.format(case_name, day_str(START_YEAR), day_str(START_MONTH))
    precip_file = nc4.Dataset(join(PRECIP_OUTPUT_DIR, file_name), 'r')
    nbins = len(precip_file.dimensions['nbins'])
    bin_lower_bounds = precip_file['bin_lower_bounds'][:]
    precip_file.close()
    ibinthresh = -1
    for i in range(nbins):
        if bin_lower_bounds[i] > PRECE_THRESHOLD:
            ibinthresh = i
            break
    if ibinthresh == -1:
        print("Warning: extreme precip threshold greater than largest bin bound.")
    prece = np.zeros((nlat, nlon))
    for imonth in imonths:
        month = months[imonth]
        year = years[imonth]
        file_name = out_file_template.format(case_name, day_str(year), day_str(month))
        precip_file = nc4.Dataset(join(PRECIP_OUTPUT_DIR, file_name), 'r')
        prece += precip_file["PRECT_amount"][:,:,ibinthresh:].sum(axis=2)
        precip_file.close()
    return prece

# Possible ways to extract a 2D section start here:
def identity(x):
    return x

def slice_at(level, x):
    return x[level,:,:]

print("Reading model output.")

varnames_readhist = [name for name in varnames if name != "PRECE"]

ref_means, test_means, diff_means = get_overall_averages(REF_CASE, TEST_CASES, months, varnames_readhist, scales)

print("Reading extreme precipitation.")
if "PRECE" in varnames:
    ref_means["PRECE"] = get_prece(REF_CASE.case_name)
    for icase in range(len(TEST_CASES)):
        test_means[icase]["PRECE"] = get_prece(TEST_CASES[icase].case_name)
        diff_means[icase]["PRECE"] = test_means[icase]["PRECE"] - ref_means["PRECE"]

# Should have this actually read from the plot_monthly_means output.
diff_global_means = {
    'CLDHGH': -0.013678454026188875,
    'CLDMED': -0.013005294090188389,
    'CLDLOW': -0.026605262766560816,
    'LWCF': -1.9367562045847697,
    'PRECC': -0.20614114936695113,
    'PRECL': 0.2874759111487673,
    'PRECT': 0.08133476183460554,
    'RELHUM': -0.6616496805116635,
    'SWCF': 5.4882057901660515,
    'TGCLDIWP': -0.34243353479029425,
    'TGCLDLWP': 0.4189776935777987,
}

for name in varnames:
    plot_name = name
    if name in plot_names:
        plot_name = plot_names[name]

    get_2D = identity
    if name in ["RELHUM", "Q", "CLDLIQ", "T", "CLOUD"]:
        get_2D = partial(slice_at, nlev-1)

    ref_plot_var = get_2D(ref_means[name])
    clim_val = [ref_plot_var.min(), ref_plot_var.max()]
    clim_diff = 0.
    for icase in range(len(TEST_CASES)):
        test_plot_var = get_2D(test_means[icase][name])
        diff_plot_var = get_2D(diff_means[icase][name])
        clim_val[0] = min(clim_val[0], test_plot_var.min())
        clim_val[1] = max(clim_val[1], test_plot_var.max())
        clim_diff = max(clim_diff, - diff_plot_var.min())
        clim_diff = max(clim_diff, diff_plot_var.max())

    if name in diff_lims:
        clim_diff = diff_lims[name]

    plt.pcolormesh(lon[:], lat[:], ref_plot_var)
    bmap.drawcoastlines()
    ax = plt.gca()
    ax.set_xticks([0., 90., 180., 270., 360.])
    ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
    ax.set_yscale('function', functions=(forward, inverse))
    ax.set_yticks([60., 45., 30., 15., 0., -15., -30., -45., -60.])
    ax.set_yticklabels(['60N', '45N', '30N', '15N', '0', '15S', '30S', '45S', '60S'])
#    ax.set_yticks([60., 30., 0., -30., -60.])
#    ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
    plt.axis('tight')
    plt.xlim([0., 360.])
    plt.colorbar()
    plt.clim(clim_val[0], clim_val[1])
    plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, REF_CASE.short_name, units[name],
                                                                  day_str(START_MONTH), day_str(START_YEAR),
                                                                  day_str(END_MONTH), day_str(END_YEAR)))
    plt.savefig('{}_{}{}.{}'.format(name, REF_CASE.short_name, suffix,
                                    PLOT_FILE_TYPE))
    plt.close()

    for icase in range(len(TEST_CASES)):
        test_plot_var = get_2D(test_means[icase][name])
        diff_plot_var = get_2D(diff_means[icase][name])
        case_name = TEST_CASES[icase].short_name

        plt.pcolormesh(lon[:], lat[:], test_plot_var)
        bmap.drawcoastlines()
        ax = plt.gca()
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yscale('function', functions=(forward, inverse))
        ax.set_yticks([60., 45., 30., 15., 0., -15., -30., -45., -60.])
        ax.set_yticklabels(['60N', '45N', '30N', '15N', '0', '15S', '30S', '45S', '60S'])
#        ax.set_yticks([60., 30., 0., -30., -60.])
#        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.axis('tight')
        plt.xlim([0., 360.])
        plt.colorbar()
        plt.clim(clim_val[0], clim_val[1])
        plt.title("{} for case {}\n({}, months {}/{} - {}/{})".format(plot_name, case_name, units[name],
                                                                      day_str(START_MONTH), day_str(START_YEAR),
                                                                      day_str(END_MONTH), day_str(END_YEAR)))
        plt.savefig('{}_{}{}.{}'.format(name, case_name, suffix,
                                        PLOT_FILE_TYPE))
        plt.close()

        plt.pcolormesh(lon[:], lat[:], diff_plot_var, cmap=cmap)
        bmap.drawcoastlines()
        ax = plt.gca()
        ax.set_xticks([0., 90., 180., 270., 360.])
        ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
        ax.set_yscale('function', functions=(forward, inverse))
        ax.set_yticks([60., 45., 30., 15., 0., -15., -30., -45., -60.])
        ax.set_yticklabels(['60N', '45N', '30N', '15N', '0', '15S', '30S', '45S', '60S'])
#        ax.set_yticks([60., 30., 0., -30., -60.])
#        ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
        plt.axis('tight')
        plt.xlim([0., 360.])
        plt.colorbar()
        plt.clim(-clim_diff, clim_diff)
        if name in diff_global_means:
            unit_string = units[name]
            if unit_string == 'fraction':
                unit_string = ''
            else:
                unit_string = " " + unit_string
            mean_string = 'mean {:.2g}'.format(diff_global_means[name]) + unit_string
        else:
            mean_string = units[name]
        plt.title("Mean difference in {} for\ncase {} ({}, months {}/{} - {}/{})".format(plot_name, case_name, mean_string,
                                                                                         day_str(START_MONTH), day_str(START_YEAR),
                                                                                         day_str(END_MONTH), day_str(END_YEAR)))
        plt.savefig('{}_diff_{}{}.{}'.format(name, case_name, suffix,
                                             PLOT_FILE_TYPE))
        plt.close()
