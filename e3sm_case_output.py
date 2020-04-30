#!/usr/bin/env python3

import numpy as np
import netCDF4 as nc4

def day_str(day):
    "Given an integer day, return the 2-digit day string used for file names."
    return "{:02d}".format(day)

class E3SMCaseOutput:

    """Class that contains information about an E3SM timestep study case.

    Methods:
    __init__
    get_daily_file_name
    get_daily_values
    compare_daily_averages

    Properties:
    case_name
    short_name
    daily_dir
    """

    def __init__(self, case_name, short_name, daily_dir):
        """Initialize a case from case directory and human-readable name.

        All arguments are used to set the corresponding properties.
        """
        self._case_name = case_name
        self._short_name = short_name
        self._daily_dir = daily_dir

    @property
    def case_name(self):
        """Name of case (as understood by CIME)."""
        return self._case_name

    @property
    def short_name(self):
        """Short human-readable name for log/plotting purposes."""
        return self._short_name

    @property
    def daily_dir(self):
        """Location of directory where daily averages are located."""
        return self._daily_dir

    def get_daily_file_name(self, day):
        """Name of the file where the average values for a given day are stored.

        Arguments:
        day - Integer representing the desired day.
        """
        return '{}/{}.0001-01-{}.nc'.format(self.daily_dir, self.case_name,
                                            day_str(day))

    def get_daily_values(self, day, varnames):
        """Retrieve the daily averages for a set of variables.

        Arguments:
        day - Integer corresponding to the desired day.
        varnames - Names of variables being read in.

        The output of this function is a dictionary associating the variable
        names to arrays containing the associated values. The time dimension is
        removed from the output arrays, but the shape is otherwise untouched.
        """
        variables = dict()
        ncfile = nc4.Dataset(self.get_daily_file_name(day), 'r')
        for varname in varnames:
            ncvar = ncfile[varname]
            var_dims = ncvar.get_dims()
            time_ind = -1
            for i in range(len(var_dims)):
                if var_dims[i].name == 'time':
                    time_ind = i
                    break
            new_dims = list(ncvar.shape)
            if time_ind != -1:
                del new_dims[time_ind]
            variables[varname] = np.zeros(new_dims)
            if time_ind == -1:
                variables[varname] = ncvar[...]
            else:
                variables[varname] = np.squeeze(ncvar, axis=time_ind)
            if "missing_value" in ncvar.ncattrs():
                miss = ncvar.missing_value
                # Following appears not to work due to rounding errors.
                # variables[varname] += np.where(variables[varname] == miss, 0., variables[varname])
                # Kludge for optical depths.
                variables[varname] = np.where(variables[varname] > 1.e35, 0., variables[varname])
        ncfile.close()
        return variables

    def compare_daily_averages(self, test_cases, day, varnames):
        """Compares the daily averages of this case to a number of other cases.

        Arguments:
        test_cases - Other case objects to compare this case to.
        day - Integer corresponding to the desired day.
        varnames - Names of variables being read in.

        The output of this function is a tuple of size 3. The first output is a
        dictionary containing the average values for this case, as in the
        get_daily_values method. The second output is a list of dictionaries
        containing the output values for each test case, in the same order as
        they were given in the test_cases argument. The third output is a list
        of dictionaries corresponding to the differences between each test case
        and this case.
        """
        ref_time_avg = self.get_daily_values(day, varnames)
        test_time_avgs = []
        diff_time_avgs = []
        for i in range(len(test_cases)):
            test_time_avgs.append(test_cases[i].get_daily_values(day, varnames))
            next_diff_time = dict()
            for varname in varnames:
                next_diff_time[varname] = test_time_avgs[i][varname] - ref_time_avg[varname]
            diff_time_avgs.append(next_diff_time)
        return (ref_time_avg, test_time_avgs, diff_time_avgs)
