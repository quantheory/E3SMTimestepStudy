#!/usr/bin/env python

from os.path import join

import numpy as np
import netCDF4 as nc4

from e3sm_case_output import day_str, time_str

CASE_NAME = "timestep_all_10s"
CASE_DIR = "/p/lscratchh/santos36/ACME/{}/run".format(CASE_NAME)
OUTPUT_DIR = "/p/lustre2/santos36/timestep_precip/"
NUM_BINS = 101
BINS_BOUNDS = (-2., 3.) # Bins between 10^-2 and 10^3 mm/day of precip.

bins = np.logspace(BINS_BOUNDS[0], BINS_BOUNDS[1], NUM_BINS-1)

START_YEAR = 4
START_MONTH = 1
END_YEAR = 4
END_MONTH = 2

month_days = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

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

filename_template = "{}.cam.h0.{}-{}-{}-{}.nc"

first_file_name = filename_template.format(CASE_NAME, "00"+day_str(START_YEAR),
                                           day_str(START_MONTH), day_str(1),
                                           time_str(0))
first_file = nc4.Dataset(join(CASE_DIR, first_file_name), 'r')
ncol = len(first_file.dimensions['ncol'])
lat = first_file['lat'][:]
lon = first_file['lon'][:]
area = first_file['area'][:]
first_file.close()

prec_vars = ("PRECC", "PRECL", "PRECSC", "PRECSL")

for i in range(nmonths):
    year = years[i]
    year_string = "00" + day_str(year)
    month = months[i]
    month_string = day_str(month)

    print("On year {}, month {}.".format(year, month))

    ndays = month_days[month-1]

    out_file_template = "{}.freq.{}-{}.nc"

    out_file_name = out_file_template.format(CASE_NAME, year_string, month_string)

    out_file = nc4.Dataset(join(OUTPUT_DIR, out_file_name), 'w')
    out_file.createDimension("ncol", ncol)
    out_file.createDimension("nbins", NUM_BINS)

    out_file.createVariable("lat", 'f8', ("ncol",))
    out_file.variables["lat"] = lat
    out_file.createVariable("lon", 'f8', ("ncol",))
    out_file.variables["lon"] = lon
    out_file.createVariable("area", 'f8', ("ncol",))
    out_file.variables["area"] = area

    out_file.createVariable("bin_lower_bounds", 'f8', ("nbins",))
    out_file.variables["bin_lower_bounds"][0] = 0.
    out_file.variables["bin_lower_bounds"][1:] = bins[:]

    var_dict = {}
    for varname in prec_vars:
        num_name = "{}_num".format(varname)
        out_file.createVariable(num_name, 'u4', ("ncol", "nbins"))
        out_file[num_name].units = "1"
        amount_name = "{}_amount".format(varname)
        out_file.createVariable(amount_name, 'f8', ("ncol", "nbins"))
        out_file[amount_name].units = "mm/day"

        var_dict[num_name] = np.zeros((ncol, NUM_BINS), dtype = np.uint32)
        var_dict[amount_name] = np.zeros((ncol, NUM_BINS))

    out_file.sample_num = ndays * 24

    for day in range(1, ndays + 1):
        print("On day {}.".format(day))
        day_string = day_str(day)
        for hour in range(0, 24):
            time_string = time_str(hour * 3600)
            in_file_name = filename_template.format(CASE_NAME, year_string,
                                                    month_string, day_string,
                                                    time_string)
            in_file = nc4.Dataset(join(CASE_DIR, in_file_name), 'r')

            for varname in prec_vars:
                var = in_file[varname][0,:]
                num_name = "{}_num".format(varname)
                amount_name = "{}_amount".format(varname)
                for i in range(ncol):
                    bin_idx = 0
                    for n in range(NUM_BINS-1):
                        if bins[n] > var[i]:
                            break
                        bin_idx += 1
                    var_dict[num_name][i, bin_idx] += 1
                    var_dict[amount_name][i, bin_idx] += var[i]

            in_file.close()

    for varname in prec_vars:
        num_name = "{}_num".format(varname)
        amount_name = "{}_amount".format(varname)
        out_file.variables[num_name] = var_dict[num_name]
        out_file.variables[amount_name] = var_dict[amount_name]

    out_file.close()
