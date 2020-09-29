#!/usr/bin/env python

from os.path import join

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

from e3sm_case_output import day_str

CASE_NAME = "timestep_all_10s"
OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip/"

TROPICS_ONLY = False

START_YEAR = 4
START_MONTH = 1
END_YEAR = 4
END_MONTH = 1

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

month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
month_weights = [month_days[months[imonth] - 1] for imonth in imonths]
weight_norm = sum(month_weights)
month_weights = [weight / weight_norm for weight in month_weights]

out_file_template = "{}.freq.{}-{}.nc"

first_file_name = out_file_template.format(CASE_NAME, "00"+day_str(START_YEAR),
                                           day_str(START_MONTH))

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

sample_num_total = 0

prec_vars = ("PRECC", "PRECL", "PRECT")

num_avgs = {}
amount_avgs = {}
for var in prec_vars:
    num_avgs[var] = np.zeros((nbins,))
    amount_avgs[var] = np.zeros((nbins,))

for i in range(nmonths):
    year = years[i]
    year_string = "00" + day_str(year)
    month = months[i]
    month_string = day_str(month)

    print("On year {}, month {}.".format(year, month))

    out_file_name = out_file_template.format(CASE_NAME, year_string, month_string)
    out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'r')

    sample_num_total += out_file.sample_num

    for var in prec_vars:
        num_name = "{}_num".format(var)
        amount_name = "{}_amount".format(var)
        for j in range(ncol):
            num_avgs[var] += out_file[num_name][j,:] * weights[j]
        for j in range(ncol):
            amount_avgs[var] += out_file[amount_name][j,:] * weights[j]

for var in prec_vars:
    num_avgs[var] /= sample_num_total
    amount_avgs[var] /= sample_num_total

for var in prec_vars:
    # Leave out zero bin from loglog plot.
    plt.loglog(bin_lower_bounds[1:], num_avgs[var][1:])
    plt.title("Frequency distribution of precipitation ({}/{}-{}/{})".format(
        day_str(START_MONTH), day_str(START_YEAR),
        day_str(END_MONTH), day_str(END_YEAR)))
    plt.xlabel("Precipitation amount (mm/day)")
    plt.ylabel("fraction")
    plt.savefig("{}_freq{}.png".format(var, suffix))
    plt.close()

    plt.semilogx(bin_lower_bounds[1:], amount_avgs[var][1:])
    plt.title("Amounts of precipitation ({}/{}-{}/{})".format(
        day_str(START_MONTH), day_str(START_YEAR),
        day_str(END_MONTH), day_str(END_YEAR)))
    plt.xlabel("Precipitation amount (mm/day)")
    plt.ylabel("Average amount (mm/day)")
    plt.savefig("{}_amount{}.png".format(var, suffix))
    plt.close()
