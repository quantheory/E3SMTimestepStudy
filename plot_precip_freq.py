#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str

REF_CASE_NAME = "timestep_ctrl"
TEST_CASE_NAME = "timestep_all_10s"
OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip/"

TROPICS_ONLY = False

START_YEAR = 3
START_MONTH = 3
END_YEAR = 4
END_MONTH = 2

suffix = '_y{}m{}-y{}m{}'.format(day_str(START_YEAR),
                                 day_str(START_MONTH),
                                 day_str(END_YEAR),
                                 day_str(END_MONTH))

if TROPICS_ONLY:
    suffix += '_tropics'

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

out_file_template = "{}.freq.{}-{}.nc"

first_file_name = out_file_template.format(REF_CASE_NAME, "00"+day_str(START_YEAR),
                                           day_str(START_MONTH))

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
    for i in range(ncol):
        if np.abs(lat[i]) > 20.:
            area[i] = 0.
area_sum = area.sum()
weights = area/area_sum
first_file.close()

ref_sample_num_total = 0
test_sample_num_total = 0

prec_vars = ("PRECC", "PRECL", "PRECT")

ref_num_avgs = {}
ref_amount_avgs = {}
for var in prec_vars:
    ref_num_avgs[var] = np.zeros((nbins,))
    ref_amount_avgs[var] = np.zeros((nbins,))

test_num_avgs = {}
test_amount_avgs = {}
for var in prec_vars:
    test_num_avgs[var] = np.zeros((nbins,))
    test_amount_avgs[var] = np.zeros((nbins,))

for i in range(nmonths):
    year = years[i]
    year_string = "00" + day_str(year)
    month = months[i]
    month_string = day_str(month)

    print("On year {}, month {}.".format(year, month))

    out_file_name = out_file_template.format(REF_CASE_NAME, year_string, month_string)
    out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'r')

    ref_sample_num_total += out_file.sample_num

    for var in prec_vars:
        num_name = "{}_num".format(var)
        amount_name = "{}_amount".format(var)
        for j in range(ncol):
            ref_num_avgs[var] += out_file[num_name][j,:] * weights[j]
        for j in range(ncol):
            ref_amount_avgs[var] += out_file[amount_name][j,:] * weights[j]

    out_file_name = out_file_template.format(TEST_CASE_NAME, year_string, month_string)
    out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'r')

    test_sample_num_total += out_file.sample_num

    for var in prec_vars:
        num_name = "{}_num".format(var)
        amount_name = "{}_amount".format(var)
        for j in range(ncol):
            test_num_avgs[var] += out_file[num_name][j,:] * weights[j]
        for j in range(ncol):
            test_amount_avgs[var] += out_file[amount_name][j,:] * weights[j]

for var in prec_vars:
    ref_num_avgs[var] /= ref_sample_num_total
    ref_amount_avgs[var] /= ref_sample_num_total
    test_num_avgs[var] /= test_sample_num_total
    test_amount_avgs[var] /= test_sample_num_total

for var in prec_vars:
    # Leave out zero bin from loglog plot.
    plt.loglog(bin_lower_bounds[1:], ref_num_avgs[var][1:], 'k')
    plt.loglog(bin_lower_bounds[1:], test_num_avgs[var][1:], 'r')
    plt.title("Frequency distribution of precipitation ({}/{}-{}/{})".format(
        day_str(START_MONTH), day_str(START_YEAR),
        day_str(END_MONTH), day_str(END_YEAR)))
    plt.xlabel("Precipitation intensity (mm/day)")
    plt.ylabel("fraction")
    plt.savefig("{}_freq{}.png".format(var, suffix))
    plt.close()

    plt.semilogx(bin_lower_bounds[1:], ref_amount_avgs[var][1:] / bin_width, 'k')
    plt.semilogx(bin_lower_bounds[1:], test_amount_avgs[var][1:] / bin_width, 'r')
    plt.title("Amounts of precipitation ({}/{}-{}/{})".format(
        day_str(START_MONTH), day_str(START_YEAR),
        day_str(END_MONTH), day_str(END_YEAR)))
    plt.xlabel("Precipitation intensity (mm/day)")
    plt.ylabel("Average precipitation amount (mm/day)")
    plt.savefig("{}_amount{}.png".format(var, suffix))
    plt.close()
