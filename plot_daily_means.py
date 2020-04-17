#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

START_DAY = 1
END_DAY = 30

DAILY_FILE_LOC="/p/lscratchh/santos36/timestep_daily_avgs/"

days = list(range(START_DAY, END_DAY+1))
ndays = len(days)

def day_str(day):
    "Given an integer day, return the 2-digit day string used for file names."
    return "{:02d}".format(day)

suffix = '_d{}-{}_presaer'.format(day_str(START_DAY), day_str(END_DAY))

log_file = open("plot_daily_log{}.txt".format(suffix), 'w')

class E3SMCaseOutput:

    """Class that contains information about an E3SM timestep study case.

    Methods:
    __init__
    get_daily_file_name
    get_daily_values

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

#REF_CASE = E3SMCaseOutput("timestep_ctrl", "CTRL", DAILY_FILE_LOC)
REF_CASE = E3SMCaseOutput("timestep_presaer_ctrl", "CTRLPA", DAILY_FILE_LOC)

TEST_CASES = [
#    E3SMCaseOutput("timestep_all_10s", "ALL10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_dyn_10s", "DYN10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_MG2_10s", "MICRO10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_CLUBB_10s", "CLUBB10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_CLUBB_MG2_10s", "CLUBBMICRO10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_CLUBB_MG2_60s", "CLUBBMICRO60", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_ZM_10s", "ZM10", DAILY_FILE_LOC),
#    E3SMCaseOutput("timestep_ZM_300s", "ZM300", DAILY_FILE_LOC),
    E3SMCaseOutput("timestep_presaer_all_10s", "ALL10PA", DAILY_FILE_LOC),
    E3SMCaseOutput("timestep_presaer_CLUBB_MG2_10s", "CLUBBMICRO10PA", DAILY_FILE_LOC),
    E3SMCaseOutput("timestep_presaer_ZM_10s", "ZM10PA", DAILY_FILE_LOC),
]

case_num = len(TEST_CASES)

rfile0 = nc4.Dataset(REF_CASE.get_daily_file_name(START_DAY), 'r')
ncol = len(rfile0.dimensions['ncol'])
area = rfile0['area'][:]
area_sum = area.sum()
weights = area/area_sum
rfile0.close()

def calc_2D_var_stats(ref_case, test_cases, day, varnames):
    ref_time_avg = ref_case.get_daily_values(day, varnames)
    test_time_avgs = []
    diff_time_avgs = []
    for i in range(len(test_cases)):
        test_time_avgs.append(test_cases[i].get_daily_values(day, varnames))
        next_diff_time = dict()
        for varname in varnames:
            next_diff_time[varname] = test_time_avgs[i][varname] - ref_time_avg[varname]
        diff_time_avgs.append(next_diff_time)
    ref_avg = dict()
    test_avgs = dict()
    diff_avgs = dict()
    rmses = dict()
    for varname in varnames:
        ref_avg[varname] = (ref_time_avg[varname] * weights).sum()
        test_avgs[varname] = []
        diff_avgs[varname] = []
        rmses[varname] = []
        for i in range(len(test_cases)):
            test_avgs[varname].append((test_time_avgs[i][varname] * weights).sum())
            diff_avgs[varname].append((diff_time_avgs[i][varname] * weights).sum())
            assert np.isclose(diff_avgs[varname][i], test_avgs[varname][i] - ref_avg[varname]), \
                "Problem with diff of variable {} from case {}".format(varname, TEST_CASE_NAMES[i])
            rmses[varname].append(np.sqrt((diff_time_avgs[i][varname]**2 * weights).sum()))
    return (ref_avg, test_avgs, diff_avgs, rmses)

def plot_vars_over_time(names, units, scales, log_plot_names):
    ref_means = dict()
    test_means = dict()
    diff_means = dict()
    rmses = dict()
    for name in names:
        ref_means[name] = np.zeros((ndays,))
        test_means[name] = np.zeros((case_num, ndays))
        diff_means[name] = np.zeros((case_num, ndays))
        rmses[name] = np.zeros((case_num, ndays))

    for iday in range(ndays):
        day = days[iday]
        print("On day: ", day, file=log_file, flush=True)
        ref_mean, test_case_means, diff_case_means, case_rmses = calc_2D_var_stats(REF_CASE, TEST_CASES, day, names)
        for name in names:
            ref_means[name][iday] = ref_mean[name]*scales[name]
            for i in range(case_num):
                test_means[name][i,iday] = test_case_means[name][i]*scales[name]
                diff_means[name][i,iday] = diff_case_means[name][i]*scales[name]
                rmses[name][i,iday] = case_rmses[name][i]*scales[name]

    for name in names:
        if name in log_plot_names:
            plot_var = plt.semilogy
        else:
            plot_var = plt.plot
        plot_var(days, ref_means[name], label=REF_CASE.short_name)
        for i in range(case_num):
            plot_var(days, test_means[name][i,:], label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} ({})".format(name, units[name]))
        plt.legend()
        plt.savefig('{}_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            plot_var(days, diff_means[name][i,:], label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("Mean {} difference ({})".format(name, units[name]))
        plt.legend()
        plt.savefig('{}_diff_time{}.png'.format(name, suffix))
        plt.close()

        for i in range(case_num):
            plt.plot(days, rmses[name][i,:], label=TEST_CASES[i].short_name)
        plt.axis('tight')
        plt.xlabel("day")
        plt.ylabel("{} RMSE ({})".format(name, units))
        plt.legend()
        plt.savefig('{}_rmse_time{}.png'.format(name, suffix))
        plt.close()

        print(name, " has reference mean: ", sum(ref_means[name])/ndays,
              file=log_file)
        for i in range(case_num):
            print(name, " has case ", TEST_CASES[i].short_name, " mean: ", sum(test_means[name][i])/ndays,
                  file=log_file)
            print(name, " has difference mean: ", sum(diff_means[name][i])/ndays,
                  file=log_file)

units = {
    'LWCF': r'$W/m^2$',
    'SWCF': r'$W/m^2$',
    'PRECC': r'$mm/day$',
    'PRECL': r'$mm/day$',
    'TGCLDCWP': r'$kg/m^2$',
    'TGCLDLWP': r'$kg/m^2$',
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
}
names = list(units.keys())
scales = dict()
for name in names:
    scales[name] = 1.
scales['PRECC'] = 1000.*86400.
scales['PRECL'] = 1000.*86400.

log_plot_names = []#'AODABS', 'AODVIS', 'AODUV']

plot_vars_over_time(names, units, scales, log_plot_names)

log_file.close()
