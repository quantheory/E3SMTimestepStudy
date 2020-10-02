#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str

REF_CASE_NAME = "timestep_all_300s"
TEST_CASE_NAMES = ["timestep_CLUBB_MG2_10s", "timestep_all_60s"]
SHORT_TEST_CASE_NAMES = ["CLUBBMICRO10", "ALL60"]
OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip/"

USE_PRESAER = False
TROPICS_ONLY = True

START_DAY = 3
END_DAY = 15

if USE_PRESAER:
    STYLES = {
        "CLUBBMICRO10PA": ('indigo', '-'),
        "ALL10PA": ('dimgrey', '-'),
        "ZM10PA": ('g', '-'),
        "CLUBBMICRO10ZM10PA": ('saddlebrown', '-'),
        "CLD10PA": ('slateblue', '-'),
    }
else:
    STYLES = {
        "DYN10": ('y', '-'),
        "CLUBB10": ('b', '-'),
        "MICRO10": ('r', '-'),
        "CLUBB10MICRO10": ('maroon', '-'),
        "CLUBBMICRO60": ('indigo', '--'),
        "CLUBBMICRO10": ('indigo', '-'),
        "ALL10": ('dimgrey', '-'),
        "ALL60": ('dimgrey', '--'),
        "ALL300": ('dimgrey', ':'),
        "ALLRAD10": ('orange', '-'),
    }

num_tests = len(TEST_CASE_NAMES)

suffix = '_d{}-d{}'.format(day_str(START_DAY), day_str(END_DAY))

if TROPICS_ONLY:
    suffix += '_tropics'

out_file_template = "{}.freq.short.d{}-d{}.nc"

first_file_name = out_file_template.format(REF_CASE_NAME, day_str(START_DAY),
                                           day_str(END_DAY))

first_file = nc4.Dataset(join(OUTPUT_DIR, first_file_name), 'r')
ncol = len(first_file.dimensions['ncol'])
nbins = len(first_file.dimensions['nbins'])
bin_lower_bounds = first_file['bin_lower_bounds'][:]
lat = first_file['lat'][:]
lon = first_file['lon'][:]
area = first_file['area'][:]
# For tropics_only cases, just use a weight of 0 for all other columns.
if TROPICS_ONLY:
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

for var in prec_vars:
    # Leave out zero bin from loglog plot.
    plt.loglog(bin_lower_bounds[1:], ref_num_avgs[var][1:], 'k')
    for i in range(num_tests):
        plt.loglog(bin_lower_bounds[1:], test_num_avgs[i][var][1:],
                   color=STYLES[SHORT_TEST_CASE_NAMES[i]][0],
                   linestyle=STYLES[SHORT_TEST_CASE_NAMES[i]][1])
    plt.title("Frequency distribution of precipitation (days {}-{})".format(
        day_str(START_DAY), day_str(END_DAY)))
    plt.xlabel("Precipitation amount (mm/day)")
    plt.ylabel("fraction")
    plt.savefig("{}_freq{}.png".format(var, suffix))
    plt.close()

    plt.semilogx(bin_lower_bounds[1:], ref_amount_avgs[var][1:], 'k')
    for i in range(num_tests):
        plt.semilogx(bin_lower_bounds[1:], test_amount_avgs[i][var][1:],
                     color=STYLES[SHORT_TEST_CASE_NAMES[i]][0],
                     linestyle=STYLES[SHORT_TEST_CASE_NAMES[i]][1])
    plt.title("Amounts of precipitation (days {}-{})".format(
        day_str(START_DAY), day_str(END_DAY)))
    plt.xlabel("Precipitation amount (mm/day)")
    plt.ylabel("Average amount (mm/day)")
    plt.savefig("{}_amount{}.png".format(var, suffix))
    plt.close()
