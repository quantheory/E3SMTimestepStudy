#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str

OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip/"

FOCUS_PRECIP = True
USE_PRESAER = False
LAND_TROPICS = False
TROPICS_ONLY = False

if LAND_TROPICS:
    TROPICS_ONLY = True

assert not (FOCUS_PRECIP and USE_PRESAER), \
    "no precipitation-specific prescribed aerosol run set has been defined"

START_DAY = 3
END_DAY = 15

if USE_PRESAER:
    REF_CASE_NAME = "timestep_presaer_ctrl"
    TEST_CASE_NAMES = [
        "timestep_presaer_ZM_10s",
        "timestep_presaer_CLUBB_MG2_10s",
        "timestep_presaer_CLUBB_MG2_10s_ZM_10s",
        "timestep_presaer_cld_10s",
        "timestep_presaer_all_10s",
        "timestep_presaer_ZM_10s_lower_tau",
        "timestep_presaer_CLUBB_MG2_10s_ZM_10s_lower_tau",
        "timestep_presaer_cld_10s_lower_tau",
        "timestep_presaer_all_10s_lower_tau",
    ]
    SHORT_TEST_CASE_NAMES = [
        "ZM10PA",
        "CLUBBMICRO10PA",
        "CLUBBMICRO10ZM10PA",
        "CLD10PA",
        "ALL10PA",
        "ZM10LTPA",
        "CLUBBMICRO10ZM10LTPA",
        "CLD10LTPA",
        "ALL10LTPA",
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
    }
elif FOCUS_PRECIP:
    REF_CASE_NAME = "timestep_ctrl"
    TEST_CASE_NAMES = [
        "timestep_MG2_10s",
#        "timestep_CLUBB_10s_MG2_10s",
#        "timestep_CLUBB_MG2_60s",
        "timestep_CLUBB_MG2_10s",
#        "timestep_all_10s",
#        "timestep_all_300s",
        "timestep_precip_grad",
        "timestep_precip_grad_MG2_10s",
        "timestep_precip_grad_CLUBB_MG2_10s",
    ]
    SHORT_TEST_CASE_NAMES = [
        "MICRO10",
#        "CLUBB10MICRO10",
#        "CLUBBMICRO60",
        "CLUBBMICRO10",
#        "ALL10",
#        "ALL300",
        "PFMG",
        "PFMGMICRO10",
        "PFMGCLUBBMICRO10",
    ]
    STYLES = {
        "MICRO10": ('r', '-'),
#        "CLUBB10MICRO10": ('maroon', '-'),
#        "CLUBBMICRO60": ('indigo', '--'),
        "CLUBBMICRO10": ('indigo', '-'),
#        "ALL10": ('dimgrey', '-'),
#        "ALL300": ('dimgrey', ':'),
        "PFMG": ('k', '-.'),
        "PFMGMICRO10": ('r', '-.'),
        "PFMGCLUBBMICRO10": ('indigo', '-.'),
    }
else:
    REF_CASE_NAME = "timestep_ctrl"
    TEST_CASE_NAMES = [
        "timestep_dyn_10s",
        "timestep_CLUBB_10s",
        "timestep_MG2_10s",
        "timestep_CLUBB_10s_MG2_10s",
        "timestep_CLUBB_MG2_Strang",
        "timestep_CLUBB_MG2_Strang_60s",
        "timestep_CLUBB_MG2_60s",
        "timestep_CLUBB_MG2_10s",
        "timestep_all_10s",
        "timestep_all_60s",
        "timestep_all_300s",
        "timestep_all_rad_10s",
    ]
    SHORT_TEST_CASE_NAMES = [
        "DYN10",
        "CLUBB10",
        "MICRO10",
        "CLUBB10MICRO10",
        "CLUBBMICROSTR",
        "CLUBBMICROSTR60",
        "CLUBBMICRO60",
        "CLUBBMICRO10",
        "ALL10",
        "ALL60",
        "ALL300",
        "ALLRAD10",
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

num_tests = len(TEST_CASE_NAMES)

suffix = '_d{}-d{}'.format(day_str(START_DAY), day_str(END_DAY))

if FOCUS_PRECIP:
    suffix += '_precip'
if USE_PRESAER:
    suffix += '_presaer'
if TROPICS_ONLY:
    if LAND_TROPICS:
        suffix += '_lndtropics'
    else:
        suffix += '_tropics'

log_file = open("plot_precip_log{}.txt".format(suffix), 'w')

out_file_template = "{}.freq.short.d{}-d{}.nc"

first_file_name = out_file_template.format(REF_CASE_NAME, day_str(START_DAY),
                                           day_str(END_DAY))

first_file = nc4.Dataset(join(OUTPUT_DIR, first_file_name), 'r')
ncol = len(first_file.dimensions['ncol'])
nbins = len(first_file.dimensions['nbins'])
bin_lower_bounds = first_file['bin_lower_bounds'][:]
bin_width = np.log(bin_lower_bounds[2] / bin_lower_bounds[1])
lat = first_file['lat'][:]
lon = first_file['lon'][:]
area = first_file['area'][:]
# For tropics_only cases, just use a weight of 0 for all other columns.
if TROPICS_ONLY:
    if LAND_TROPICS:
        # Just pick a random file with the same grid as the run.
        landfrac_file_name = '/p/lustre2/santos36/timestep_monthly_avgs/timestep_ctrl.0001-01.nc'
        landfrac_file = nc4.Dataset(landfrac_file_name, 'r')
        landfrac = landfrac_file['LANDFRAC'][0,:]
        for i in range(ncol):
            if np.abs(lat[i]) > 30.:
                area[i] = 0.
            else:
                area[i] *= landfrac[i]
        landfrac_file.close()
    else:
        for i in range(ncol):
            if np.abs(lat[i]) > 30.:
                area[i] = 0.
area_sum = area.sum()
weights = area/area_sum
first_file.close()

ref_sample_num_total = 0
test_sample_num_totals = [0 for i in range(num_tests)]

prec_vars = ("PRECC", "PRECL", "PRECT")

ref_num_avgs = {}
ref_amount_avgs = {}
for var in prec_vars:
    ref_num_avgs[var] = np.zeros((nbins,))
    ref_amount_avgs[var] = np.zeros((nbins,))

test_num_avgs = [{} for i in range (num_tests)]
test_amount_avgs = [{} for i in range (num_tests)]
for i in range(num_tests):
    for var in prec_vars:
        test_num_avgs[i][var] = np.zeros((nbins,))
        test_amount_avgs[i][var] = np.zeros((nbins,))

out_file_name = out_file_template.format(REF_CASE_NAME, day_str(START_DAY),
                                         day_str(END_DAY))
out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'r')

ref_sample_num_total += out_file.sample_num

for var in prec_vars:
    num_name = "{}_num".format(var)
    amount_name = "{}_amount".format(var)
    for j in range(ncol):
        ref_num_avgs[var] += out_file[num_name][j,:] * weights[j]
    for j in range(ncol):
        ref_amount_avgs[var] += out_file[amount_name][j,:] * weights[j]

for i in range(num_tests):
    out_file_name = out_file_template.format(TEST_CASE_NAMES[i], day_str(START_DAY),
                                             day_str(END_DAY))
    out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'r')

    test_sample_num_totals[i] += out_file.sample_num

    for var in prec_vars:
        num_name = "{}_num".format(var)
        amount_name = "{}_amount".format(var)
        for j in range(ncol):
            test_num_avgs[i][var] += out_file[num_name][j,:] * weights[j]
        for j in range(ncol):
            test_amount_avgs[i][var] += out_file[amount_name][j,:] * weights[j]

for var in prec_vars:
    ref_num_avgs[var] /= ref_sample_num_total
    ref_amount_avgs[var] /= ref_sample_num_total
    for i in range(num_tests):
        test_num_avgs[i][var] /= test_sample_num_totals[i]
        test_amount_avgs[i][var] /= test_sample_num_totals[i]

# Threshold for precipitation to be considered "extreme", in mm/day.
PRECE_THRESHOLD = 97.
ibinthresh = -1
for i in range(nbins):
    if bin_lower_bounds[i] > PRECE_THRESHOLD:
        ibinthresh = i
        break
if ibinthresh == -1:
    print("Warning: extreme precip threshold greater than largest bin bound.")

for var in prec_vars:
    # Leave out zero bin from loglog plot.
    plt.loglog(bin_lower_bounds[1:], ref_num_avgs[var][1:], 'k')
    for i in range(num_tests):
        plt.loglog(bin_lower_bounds[1:], test_num_avgs[i][var][1:],
                   color=STYLES[SHORT_TEST_CASE_NAMES[i]][0],
                   linestyle=STYLES[SHORT_TEST_CASE_NAMES[i]][1])
    plt.title("Frequency distribution of precipitation (days {}-{})".format(
        day_str(START_DAY), day_str(END_DAY)))
    plt.xlabel("Precipitation intensity (mm/day)")
    plt.ylabel("fraction")
    plt.savefig("{}_freq{}.png".format(var, suffix))
    plt.close()

    plt.semilogx(bin_lower_bounds[1:], ref_amount_avgs[var][1:] / bin_width, 'k')
    if var == "PRECT":
        print("Extreme precipitation rate for reference: ",
              ref_amount_avgs[var][ibinthresh:].sum(),
              file=log_file)
    for i in range(num_tests):
        plt.semilogx(bin_lower_bounds[1:], test_amount_avgs[i][var][1:] / bin_width,
                     color=STYLES[SHORT_TEST_CASE_NAMES[i]][0],
                     linestyle=STYLES[SHORT_TEST_CASE_NAMES[i]][1])
        if var == "PRECT":
            print("Extreme precipitation rate for ", SHORT_TEST_CASE_NAMES[i], ": ",
                  test_amount_avgs[i][var][ibinthresh:].sum(), "(Diff = ",
                  test_amount_avgs[i][var][ibinthresh:].sum() - ref_amount_avgs[var][ibinthresh:].sum(), ")",
                  file=log_file)
    plt.title("Amounts of precipitation (days {}-{})".format(
        day_str(START_DAY), day_str(END_DAY)))
    plt.xlabel("Precipitation intensity (mm/day)")
    plt.ylabel("Average precipitation amount (mm/day)")
    plt.savefig("{}_amount{}.png".format(var, suffix))
    plt.close()

log_file.close()
