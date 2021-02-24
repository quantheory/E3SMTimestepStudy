#!/usr/bin/env python3

from functools import partial

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import E3SMCaseOutput, day_str

START_DAY = 3
END_DAY = 30

DAILY_FILE_LOC="/p/lscratchh/santos36/timestep_daily_avgs/"

USE_PRESAER=False
LAND_TROPICS = False
TROPICS_ONLY=True
MIDLATITUDES_ONLY=False

if LAND_TROPICS:
    TROPICS_ONLY = True

assert not (TROPICS_ONLY and MIDLATITUDES_ONLY), \
    "can't do only tropics and only midlatitudes"

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)

suffix = '_d{}-{}'.format(day_str(START_DAY), day_str(END_DAY))

if USE_PRESAER:
    suffix += '_presaer'
if TROPICS_ONLY:
    if LAND_TROPICS:
        suffix += '_lndtropics'
    else:
        suffix += '_tropics'
if MIDLATITUDES_ONLY:
    suffix += '_midlats'

log_file = open("plot_dynamics_attribution_daily_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
    ]
    STYLES = {
        "CLUBBMICRO10PA": ('indigo', '-'),
        "ALL10PA": ('dimgrey', '-'),
        "ZM10PA": ('g', '-'),
        "CLUBBMICRO10ZM10PA": ('saddlebrown', '-'),
        "CLD10PA": ('slateblue', '-'),
        "ALL10LTPA": ('dimgrey', '-.'),
        "ZM10LTPA": ('g', '-.'),
        "CLUBBMICRO10ZM10LTPA": ('saddlebrown', '-.'),
        "CLD10LTPA": ('slateblue', '-.'),
        "CLD10LT2PA": ('slateblue', ':'),
    }
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_all_10s", "ALL10", DAILY_FILE_LOC, START_DAY, END_DAY),
    ]
    STYLES = {
        "DYN10": ('y', '-'),
        "CLUBB10": ('b', '-'),
        "MICRO10": ('r', '-'),
        "CLUBB10MICRO10": ('maroon', '-'),
        "CLUBBMICROSTR": ('m', '-'),
        "CLUBBMICROSTR60": ('m', '--'),
        "CLUBBMICRO60": ('indigo', '--'),
        "CLUBBMICRO10": ('indigo', '-'),
        "ALL10": ('r', '-'),
        "ALL60": ('dimgrey', '--'),
        "ALL300": ('dimgrey', ':'),
        "ALLRAD10": ('orange', '-'),
    }

case_num = len(TEST_CASES)

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
nlev = len(rfile0.dimensions['lev'])
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
# For tropics_only cases, just use a weight of 0 for all other cases.
if TROPICS_ONLY:
    lat = rfile0['lat'][:]
    if LAND_TROPICS:
        landfrac = rfile0['LANDFRAC'][0,:]
        for i in range(ncol):
            if np.abs(lat[i]) > 30.:
                area[i] = 0.
            else:
                area[i] *= landfrac[i]
    else:
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

nbins = 100
omega_span = (-2.7, 1.)
bin_size = (omega_span[1] - omega_span[0]) / nbins

def time_avg_to_bins(varnames, time_avg):
    omega_pdf = np.zeros((nbins,))
    var_dists = dict()
    for varname in varnames:
        if varname in vars_3D:
            var_dists[varname] = np.zeros((nlev, nbins))
        else:
            var_dists[varname] = np.zeros((nbins,))
    for i in range(ncol):
        omega = time_avg["OMEGA500"][i]
        ibin = int((omega - omega_span[0]) / bin_size)
        if 0 <= ibin < nbins:
            weight = weights[i]
            omega_pdf[ibin] += weight
            for varname in varnames:
                if varname in vars_3D:
                    for jlev in range(nlev):
                        var_dists[varname][jlev,ibin] += time_avg[varname][jlev,i] * weights[i]
                else:
                    var_dists[varname][ibin] += time_avg[varname][i] * weights[i]
    for varname in varnames:
        if varname in vars_3D:
            for jlev in range(nlev):
                for ibin in range(nbins):
                    if omega_pdf[ibin] > 0:
                        var_dists[varname][jlev,ibin] /= omega_pdf[ibin]
        else:
            for ibin in range(nbins):
                if omega_pdf[ibin] > 0:
                    var_dists[varname][ibin] /= omega_pdf[ibin]
    return (omega_pdf, var_dists)

def calc_omega_var_bins(ref_case, test_cases, day, varnames):
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
    ref_time_avg, test_time_avgs, diff_time_avgs = ref_case.compare_daily_averages(test_cases, day, varnames_read)
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
    ref_omega_pdf, ref_var_dists = time_avg_to_bins(varnames, ref_time_avg)
    test_omega_pdfs = []
    test_var_dists = []
    for i in range(len(test_cases)):
        test_pdf, test_dist = time_avg_to_bins(varnames, test_time_avgs[i])
        test_omega_pdfs.append(test_pdf)
        test_var_dists.append(test_dist)
    return (ref_omega_pdf, ref_var_dists, test_omega_pdfs, test_var_dists)

# Possible ways to extract a 2D section start here:
def identity(x):
    return x

def slice_at(level, x):
    return x[level,:]

def plot_dyn_attr(names, units, scales, log_plot_names):

    ref_total_pdf = np.zeros((nbins,))
    ref_total_dists = dict()
    for name in names:
        if name in vars_3D:
            ref_total_dists[name] = np.zeros((nlev, nbins))
        else:
            ref_total_dists[name] = np.zeros((nbins,))
    test_total_pdfs = []
    test_total_dists = []
    for i in range(case_num):
        test_total_pdfs.append(np.zeros((nbins,)))
        total_dists = dict()
        for name in names:
            if name in vars_3D:
                total_dists[name] = np.zeros((nlev, nbins))
            else:
                total_dists[name] = np.zeros((nbins,))
        test_total_dists.append(total_dists)

    for iday in range(ndays):
        day = days[iday]
        print("On day: ", day, file=log_file, flush=True)
        ref_omega_pdf, ref_var_dists, test_omega_pdfs, test_var_dists = calc_omega_var_bins(REF_CASE, TEST_CASES, day, names)
        ref_total_pdf += ref_omega_pdf / ndays
        for name in names:
            ref_total_dists[name] += ref_var_dists[name]*scales[name] / ndays
        for i in range(case_num):
            test_total_pdfs[i] += test_omega_pdfs[i] / ndays
            for name in names:
                test_total_dists[i][name] += test_var_dists[i][name]*scales[name] / ndays

    omegas = 864.*np.linspace(omega_span[0] + 0.5*bin_size, omega_span[1] - 0.5*bin_size, nbins)
    plt.semilogy(omegas, ref_total_pdf / bin_size, label=REF_CASE.short_name,
                 color='k')
    for i in range(case_num):
        plt.semilogy(omegas, test_total_pdfs[i] / bin_size,
                     label=TEST_CASES[i].short_name,
                     color=STYLES[TEST_CASES[i].short_name][0],
                     linestyle=STYLES[TEST_CASES[i].short_name][1])
    plt.axis('tight')
    plt.xlabel(r'$\omega_{500}$')
    plt.ylabel(r'$p(\omega_{500})$')
    plt.ylim(1.e-6, 1.e1)
    plt.savefig('OMEGA500_dist{}.png'.format(suffix))
    plt.close()

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
        ref_plot_var = get_2D(ref_total_dists[name]) * ref_total_pdf
        plot_var(omegas, ref_plot_var, label=REF_CASE.short_name, color='k')
        for i in range(case_num):
            test_plot_var = get_2D(test_total_dists[i][name]) * test_total_pdfs[i]
            plot_var(omegas,
                     test_plot_var,
                     label=TEST_CASES[i].short_name,
                     color=STYLES[TEST_CASES[i].short_name][0],
                     linestyle=STYLES[TEST_CASES[i].short_name][1])
            print(name, " for case ", TEST_CASES[i].short_name, " has total diff ",
                  (test_plot_var - ref_plot_var).sum(),
                  file=log_file, flush=True)
        plt.axis('tight')
        plt.xlabel(r'$\omega_{500}$')
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_pdf{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            pdf_diff = test_total_pdfs[i] - ref_total_pdf
            test_plot_var = get_2D(ref_total_dists[name]) * pdf_diff
            plot_var(omegas,
                     test_plot_var,
                     label=TEST_CASES[i].short_name)
            print(name, " for case ", TEST_CASES[i].short_name, " has dynamics diff ",
                  test_plot_var.sum(), file=log_file, flush=True)
        plt.axis('tight')
        plt.xlabel(r'$\omega_{500}$')
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_pdf_dyn{}.png'.format(name, suffix))
        plt.close()

        ref_plot_var = get_2D(ref_total_dists[name])
        for i in range(case_num):
            test_plot_var = (get_2D(test_total_dists[i][name]) - ref_plot_var) * ref_total_pdf
            plot_var(omegas,
                     test_plot_var,
                     label=TEST_CASES[i].short_name)
            print(name, " for case ", TEST_CASES[i].short_name, " has thermodynamics diff ",
                  test_plot_var.sum(), file=log_file, flush=True)
        plt.axis('tight')
        plt.xlabel(r'$\omega_{500}$')
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_pdf_thermo{}.png'.format(name, suffix))
        plt.close()

        ref_plot_var = get_2D(ref_total_dists[name])
        for i in range(case_num):
            pdf_diff = test_total_pdfs[i] - ref_total_pdf
            test_plot_var = (get_2D(test_total_dists[i][name]) - ref_plot_var) * pdf_diff
            plot_var(omegas,
                     test_plot_var,
                     label=TEST_CASES[i].short_name)
            print(name, " for case ", TEST_CASES[i].short_name, " has covariance diff ",
                  test_plot_var.sum(), file=log_file, flush=True)
        plt.axis('tight')
        plt.xlabel(r'$\omega_{500}$')
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_pdf_covar{}.png'.format(name, suffix))
        plt.close()

plot_names = {
    'LWCF': "long-wave cloud forcing",
    'SWCF': "short wave cloud forcing",
    'PRECC': "convective precipitation",
    'PRECL': "large scale precipitation",
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
    'TMQ': r'$kg/m^2$',
    'T': r'$K$',
}
#names = list(units.keys())
names = ['PRECC', 'PRECL', 'PRECT', 'OMEGA500']
scales = dict()
for name in names:
    scales[name] = 1.
scales['TGCLDIWP'] = 1000.
scales['TGCLDLWP'] = 1000.
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

plot_dyn_attr(names, units, scales, log_plot_names)

log_file.close()
