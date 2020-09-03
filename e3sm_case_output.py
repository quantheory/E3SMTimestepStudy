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
    day_is_available
    compare_daily_averages

    Properties:
    case_name
    short_name
    daily_dir
    """

    def __init__(self, case_name, short_name, daily_dir, start_day, end_day):
        """Initialize a case from case directory and human-readable name.

        All arguments are used to set the corresponding properties.
        """
        self._case_name = case_name
        self._short_name = short_name
        self._daily_dir = daily_dir
        self._start_day = start_day
        self._end_day = end_day

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

    @property
    def start_day(self):
        """First day for which output is available."""
        return self._start_day

    @property
    def end_day(self):
        """Last day for which output is available."""
        return self._end_day

    def day_is_available(self, day):
        """Returns True if the given day has available output for this case."""
        return self.start_day <= day and self.end_day >= day

    def get_daily_file_name(self, day):
        """Name of the file where the average values for a given day are stored.

        Arguments:
        day - Integer representing the desired day.
        """
        assert self.day_is_available(day)
        return '{}/{}.0001-01-{}.nc'.format(self.daily_dir, self.case_name,
                                            day_str(day))

    def get_monthly_file_name(self, month, year):
        """Name of the file where the average values for a given month are stored.

        Arguments:
        month - Integer representing the desired month.
        year - Integer representing year containing the desired month.
        """
        #assert self.day_is_available(day)
        return '{}/{}.00{}-{}.nc'.format(self.daily_dir, self.case_name,
                                         day_str(year), day_str(month))

    def get_daily_values(self, day, varnames):
        """Retrieve the daily averages for a set of variables.

        Arguments:
        day - Integer corresponding to the desired day.
        varnames - Names of variables being read in.

        The output of this function is a dictionary associating the variable
        names to arrays containing the associated values. The time dimension is
        removed from the output arrays, but the shape is otherwise untouched.
        """
        assert self.day_is_available(day)
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

    def get_monthly_values(self, month, year, varnames):
        """Retrieve the monthly averages for a set of variables.

        Arguments:
        month - Integer corresponding to the desired month.
        year - Integer representing year containing the desired month.
        varnames - Names of variables being read in.

        The output of this function is a dictionary associating the variable
        names to arrays containing the associated values. The time dimension is
        removed from the output arrays, but the shape is otherwise untouched.
        """
        # assert self.day_is_available(day)
        variables = dict()
        ncfile = nc4.Dataset(self.get_monthly_file_name(month, year), 'r')
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

        For test cases that do not extend to the given day, the corresponding
        item in the second and third output lists will be None.

        """
        ref_time_avg = self.get_daily_values(day, varnames)
        test_time_avgs = []
        diff_time_avgs = []
        for i in range(len(test_cases)):
            if test_cases[i].day_is_available(day):
                test_time_avgs.append(test_cases[i].get_daily_values(day, varnames))
                next_diff_time = dict()
                for varname in varnames:
                    next_diff_time[varname] = test_time_avgs[i][varname] - ref_time_avg[varname]
                diff_time_avgs.append(next_diff_time)
            else:
                test_time_avgs.append(None)
                diff_time_avgs.append(None)
        return (ref_time_avg, test_time_avgs, diff_time_avgs)

    def compare_monthly_averages(self, test_cases, month, year, varnames):
        """Compares the monthly averages of this case to a number of other cases.

        Arguments:
        test_cases - Other case objects to compare this case to.
        month - Integer corresponding to the desired month.
        year - Integer representing year containing the desired month.
        varnames - Names of variables being read in.

        The output of this function is a tuple of size 3. The first output is a
        dictionary containing the average values for this case, as in the
        get_monthly_values method. The second output is a list of dictionaries
        containing the output values for each test case, in the same order as
        they were given in the test_cases argument. The third output is a list
        of dictionaries corresponding to the differences between each test case
        and this case.

        """
        ref_time_avg = self.get_monthly_values(month, year, varnames)
        test_time_avgs = []
        diff_time_avgs = []
        for i in range(len(test_cases)):
            test_time_avgs.append(test_cases[i].get_monthly_values(month, year, varnames))
            next_diff_time = dict()
            for varname in varnames:
                next_diff_time[varname] = test_time_avgs[i][varname] - ref_time_avg[varname]
            diff_time_avgs.append(next_diff_time)
        return (ref_time_avg, test_time_avgs, diff_time_avgs)
