#!/usr/bin/env python3

from functools import partial

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

FOCUS_PRECIP = False
USE_PRESAER = True
LAND_ONLY = False
OCEAN_ONLY = False # Note that this includes ocean and sea-ice grid cells
TROPICS_ONLY = False
MIDLATITUDES_ONLY = False

assert not (FOCUS_PRECIP and USE_PRESAER), \
    "no precipitation-specific prescribed aerosol run set has been defined"

assert not (LAND_ONLY and OCEAN_ONLY), \
    "can't do only land and only ocean"

assert not (TROPICS_ONLY and MIDLATITUDES_ONLY), \
    "can't do only tropics and only midlatitudes"

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)
navgdays = END_AVG_DAY - START_AVG_DAY + 1

suffix = '_d{}-{}'.format(day_str(START_DAY), day_str(END_DAY))
if FOCUS_PRECIP:
    suffix += '_precip'
if USE_PRESAER:
    suffix += '_presaer'
if LAND_ONLY:
    sfc_suffix = 'lnd'
elif OCEAN_ONLY:
    sfc_suffix = 'ocn'
else:
    sfc_suffix = ''
if TROPICS_ONLY:
    suffix += '_{}tropics'.format(sfc_suffix)
elif MIDLATITUDES_ONLY:
    suffix += '_{}midlats'.format(sfc_suffix)
elif sfc_suffix != '':
    suffix += '_{}'.format(sfc_suffix)

log_file = open("plot_daily_log{}.txt".format(suffix), 'w')

if USE_PRESAER:
    REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_presaer_ZM_10s", "ZM10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_ZM_10s_lower_tau", "ZM10LTPA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s", "CLUBBMICRO10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s_ZM_10s", "CLUBBMICRO10ZM10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s_ZM_10s_lower_tau", "CLUBBMICRO10ZM10LTPA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_cld_10s", "CLD10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_cld_10s_lower_tau", "CLD10LTPA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_cld_10s_lower_tau2", "CLD10LT2PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_all_10s", "ALL10PA", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_presaer_all_10s_lower_tau", "ALL10LTPA", DAILY_FILE_LOC, START_DAY, END_DAY),
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
elif FOCUS_PRECIP:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_precip_grad", "PFMG", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_CLUBB_10s", "CLUBB10", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_CLUBB_10s_MG2_10s", "CLUBB10MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_MG2_10s", "MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_precip_grad_MG2_10s", "PFMGMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_precip_grad_CLUBB_MG2_10s", "PFMGCLUBBMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_all_300s", "ALL300", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_all_60s", "ALL60", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_all_10s", "ALL10", DAILY_FILE_LOC, START_DAY, END_DAY),
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
        "ALL10": ('dimgrey', '-'),
        "ALL60": ('dimgrey', '--'),
        "ALL300": ('dimgrey', ':'),
        "ALLRAD10": ('orange', '-'),
        "PFMG": ('k', '-.'),
        "PFMGMICRO10": ('r', '-.'),
        "PFMGCLUBBMICRO10": ('indigo', '-.'),
    }
else:
    REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC, START_DAY, END_DAY)
    TEST_CASES = [
        E3SMCaseOutput("timestep_dyn_10s", "DYN10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_10s", "CLUBB10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_10s_MG2_10s", "CLUBB10MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_CLUBB_MG2_Strang", "CLUBBMICROSTR", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_CLUBB_MG2_Strang_60s", "CLUBBMICROSTR60", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_MG2_10s", "MICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC, START_DAY, END_DAY),
#        E3SMCaseOutput("timestep_ZM_10s", "ZM10", DAILY_FILE_LOC, START_DAY, END_ZM10S_DAY),
#        E3SMCaseOutput("timestep_ZM_300s", "ZM300", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_rad_10s", "ALLRAD10", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_300s", "ALL300", DAILY_FILE_LOC, START_DAY, END_DAY),
        E3SMCaseOutput("timestep_all_60s", "ALL60", DAILY_FILE_LOC, START_DAY, END_DAY),
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
        "ALL10": ('dimgrey', '-'),
        "ALL60": ('dimgrey', '--'),
        "ALL300": ('dimgrey', ':'),
        "ALLRAD10": ('orange', '-'),
    }

case_num = len(TEST_CASES)

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
nlev = len(rfile0.dimensions['lev'])
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
if LAND_ONLY:
    landfrac = rfile0['LANDFRAC'][0,:]
    area *= landfrac
elif OCEAN_ONLY:
    landfrac = rfile0['LANDFRAC'][0,:]
    area *= 1. - landfrac
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

def calc_var_stats(ref_case, test_cases, day, varnames):
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
            if test_cases[i].day_is_available(day):
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
            else:
                test_avgs[varname].append(None)
                diff_avgs[varname].append(None)
                rmses[varname].append(None)
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
            ref_means[name] = np.zeros((ndays, nlev))
            test_means[name] = np.zeros((case_num, ndays, nlev))
            diff_means[name] = np.zeros((case_num, ndays, nlev))
            rmses[name] = np.zeros((case_num, ndays, nlev))
        else:
            ref_means[name] = np.zeros((ndays,))
            test_means[name] = np.zeros((case_num, ndays))
            diff_means[name] = np.zeros((case_num, ndays))
            rmses[name] = np.zeros((case_num, ndays))

    for iday in range(ndays):
        day = days[iday]
        print("On day: ", day, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_var_stats(REF_CASE, TEST_CASES, day, names)
        for name in names:
            ref_means[name][iday] = ref_mean[name]*scales[name]
            for i in range(case_num):
                if TEST_CASES[i].day_is_available(day):
                    test_means[name][i,iday] = test_case_means[name][i]*scales[name]
                    diff_means[name][i,iday] = diff_case_means[name][i]*scales[name]
                    rmses[name][i,iday] = case_rmses[name][i]*scales[name]

    for name in names:
        plot_name = name
        if name in plot_names:
            plot_name = plot_names[name]

        get_2D = identity
        if name in vars_3D:
            get_2D = partial(slice_at, nlev-1)

        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        for i in range(case_num):
            test_plot_var = get_2D(test_means[name][i])
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     test_plot_var[start_ind:end_ind],
                     label=TEST_CASES[i].short_name,
                     color=STYLES[TEST_CASES[i].short_name][0],
                     linestyle=STYLES[TEST_CASES[i].short_name][1])
        ref_plot_var = get_2D(ref_means[name])
        plot_var(days, ref_plot_var, label=REF_CASE.short_name, color='k')
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} ({})".format(plot_name, units[name]))
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            diff_plot_var = get_2D(diff_means[name][i])
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     diff_plot_var[start_ind:end_ind],
                     label=TEST_CASES[i].short_name,
                     color=STYLES[TEST_CASES[i].short_name][0],
                     linestyle=STYLES[TEST_CASES[i].short_name][1])
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} difference ({})".format(plot_name, units[name]))
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            rmse_plot_var = get_2D(rmses[name][i])
            start_ind = TEST_CASES[i].start_day - START_DAY
            end_ind = TEST_CASES[i].end_day - START_DAY + 1
            plot_var(days[start_ind:end_ind],
                     rmse_plot_var[start_ind:end_ind],
                     label=TEST_CASES[i].short_name,
                     color=STYLES[TEST_CASES[i].short_name][0],
                     linestyle=STYLES[TEST_CASES[i].short_name][1])
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("{} RMSE ({})".format(plot_name, units[name]))
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        print(name, " has reference mean: ", sum(ref_plot_var[START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
              file=log_file)
        for i in range(case_num):
            case_name = TEST_CASES[i].short_name
            test_plot_var = get_2D(test_means[name][i])
            diff_plot_var = get_2D(diff_means[name][i])
            print(name, " has case ", case_name, " mean: ", sum(test_plot_var[START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
                  file=log_file)
            print(name, " has difference mean: ", sum(diff_plot_var[START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
                  file=log_file)
            if USE_PRESAER and "LT" in case_name:
                compare_name = TEST_CASES[i-1].short_name
                compare_plot_var = get_2D(test_means[name][i-1])
                print(name, " has mean difference from ", compare_name, ": ",
                      sum(test_plot_var[START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays - \
                      sum(compare_plot_var[START_AVG_DAY-START_DAY:END_AVG_DAY-START_DAY+1])/navgdays,
                      file=log_file)

plot_names = {
    'LWCF': "longwave cloud forcing",
    'SWCF': "shortwave cloud forcing",
    'PRECC': "convective precipitation",
    'PRECL': "large-scale precipitation",
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
    'TMQ': "precipitable water",
    'CLOUD': "lowest level cloud fraction",
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
    'CLDLIQ': r"$g/kg$",
    'TMQ': r'$kg/m^2$',
    'CLOUD': r'$fraction$',
    'T': r'$K$',
}
names = list(units.keys())
scales = dict()
for name in names:
    scales[name] = 1.
scales['TGCLDIWP'] = 1000.
scales['TGCLDLWP'] = 1000.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.
scales['PRECT'] = 1000.*86400.
scales['Q'] = 1000.
scales['CLDLIQ'] = 1000.

vars_3D = [
    'RELHUM',
    'Q',
    'CLDLIQ',
    'T',
    'CLOUD',
]

log_plot_names = []#'AODABS', 'AODVIS', 'AODUV']

plot_vars_over_time(names, units, scales, log_plot_names)

log_file.close()
