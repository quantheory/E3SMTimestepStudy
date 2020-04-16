#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import netCDF4 as nc4

CONTROL_FILE_NAME = "/home/santos/Data/GATEIII_ctrl.cam.h0.1974-08-30-00000.nc"
TEST1_FILE_NAME = "/home/santos/Data/GATEIII_ZM_10s.cam.h0.1974-08-30-00000.nc"
TEST2_FILE_NAME = "/home/santos/Data/GATEIII_CLUBB_MG2_10s.cam.h0.1974-08-30-00000.nc"
CONTROL_DESC = "Control"
TEST1_DESC = "ZM at 10s"
TEST2_DESC = "CLUBB and MG2 at 10s"
SUFFIX_CTRL = "_ctrl"
SUFFIX1="_ZM_10s"
SUFFIX2="_CLUBB_MG2_10s"
PLOT_VARS = ['SWCF', 'PRECC', 'PRECL']
STEPS_IN_AVG = {
    'SWCF': 48,
    'PRECC': 6,
    'PRECL': 6,
}
PLOT_VERT_VARS = ['CLDLIQ']
STEPS_IN_VERT_AVG = {
    'CLDLIQ': 1,
}
DIFF_CMAP = plt.get_cmap('coolwarm')

# GATEIII specific variables
NUM_DAYS = 20
STEPS_PER_DAY = 48

total_steps = NUM_DAYS*STEPS_PER_DAY

ctrl_file = nc4.Dataset(CONTROL_FILE_NAME, 'r')
test1_file = nc4.Dataset(TEST1_FILE_NAME, 'r')
test2_file = nc4.Dataset(TEST2_FILE_NAME, 'r')

times = ctrl_file.variables['time'][:]
time_len = len(times)

nlev = len(ctrl_file.dimensions['lev'])
lev = ctrl_file.variables['lev'][:]

for var_name, avg_step in zip(PLOT_VARS, STEPS_IN_AVG):
    assert total_steps % avg_step == 0
    num_avg = total_steps // avg_step
    plot_times = np.linspace(0, NUM_DAYS+1, num_avg+1)

    var_units = ctrl_file.variables[var_name].units

    ctrl_var = ctrl_file.variables[var_name][:,0]
    test1_var = test1_file.variables[var_name][:,0]
    test2_var = test2_file.variables[var_name][:,0]

    filename = "{}_compare.png".format(var_name)

    ctrl_avg = np.zeros((num_avg,))
    test1_avg = np.zeros((num_avg,))
    test2_avg = np.zeros((num_avg,))
    for i in range(num_avg):
        ctrl_avg[i] = ctrl_var[i*avg_step:(i+1)*avg_step].sum()
        test1_avg[i] = test1_var[i*avg_step:(i+1)*avg_step].sum()
        test2_avg[i] = test2_var[i*avg_step:(i+1)*avg_step].sum()

    plt.plot(plot_times[1:], ctrl_avg, 'k', label=CONTROL_DESC)
    plt.plot(plot_times[1:], test1_avg, 'r', label=TEST1_DESC)
    plt.plot(plot_times[1:], test2_avg, 'b', label=TEST2_DESC)
    plt.title("{} comparison".format(var_name))
    plt.xlabel("Days since model start")
    plt.ylabel("{} ({})".format(var_name, var_units))
    plt.legend(loc='best')
    plt.savefig(filename)
    plt.close()

for var_name, avg_step in zip(PLOT_VERT_VARS, STEPS_IN_VERT_AVG):
    assert total_steps % avg_step == 0
    num_avg = total_steps // avg_step
    plot_times = np.linspace(0, NUM_DAYS+1, num_avg+1)

    var_units = ctrl_file.variables[var_name].units

    ctrl_var = ctrl_file.variables[var_name][:,:,0]
    test1_var = test1_file.variables[var_name][:,:,0]
    test2_var = test2_file.variables[var_name][:,:,0]

    ctrl_filename = "{}{}.png".format(var_name, SUFFIX_CTRL)
    test1_filename = "{}{}.png".format(var_name, SUFFIX1)
    test2_filename = "{}{}.png".format(var_name, SUFFIX2)
    diff1_filename = "{}_diff{}.png".format(var_name, SUFFIX1)
    diff2_filename = "{}_diff{}.png".format(var_name, SUFFIX2)

    ctrl_avg = np.zeros((num_avg, nlev))
    test1_avg = np.zeros((num_avg, nlev))
    test2_avg = np.zeros((num_avg, nlev))
    for i in range(num_avg):
        for j in range(nlev):
            ctrl_avg[i,j] = ctrl_var[i*avg_step:(i+1)*avg_step,j].sum()
            test1_avg[i,j] = test1_var[i*avg_step:(i+1)*avg_step,j].sum()
            test2_avg[i,j] = test2_var[i*avg_step:(i+1)*avg_step,j].sum()

    diff1_avg = test1_avg - ctrl_avg
    diff2_avg = test2_avg - ctrl_avg

    plt.pcolor(plot_times, lev, ctrl_avg.T)
    plt.title("{} results from {}".format(var_name, CONTROL_DESC))
    plt.xlabel("Days since model start")
    plt.ylabel("Pressure (hPa)".format(var_name, var_units))
    plt.colorbar()
    plt.savefig(ctrl_filename)
    plt.close()

    plt.pcolor(plot_times, lev, test1_avg.T)
    plt.title("{} results from {}".format(var_name, TEST1_DESC))
    plt.xlabel("Days since model start")
    plt.ylabel("Pressure (hPa)".format(var_name, var_units))
    plt.colorbar()
    plt.savefig(test1_filename)
    plt.close()

    plt.pcolor(plot_times, lev, diff1_avg.T, cmap=DIFF_CMAP)
    plt.title("{} diff from {} to {}".format(var_name, CONTROL_DESC, TEST1_DESC))
    plt.xlabel("Days since model start")
    plt.ylabel("Pressure (hPa)".format(var_name, var_units))
    plt.colorbar()
    plt.savefig(diff1_filename)
    plt.close()

    plt.pcolor(plot_times, lev, test2_avg.T)
    plt.title("{} results from {}".format(var_name, TEST2_DESC))
    plt.xlabel("Days since model start")
    plt.ylabel("Pressure (hPa)".format(var_name, var_units))
    plt.colorbar()
    plt.savefig(test2_filename)
    plt.close()

    plt.pcolor(plot_times, lev, diff2_avg.T, cmap=DIFF_CMAP)
    plt.title("{} diff from {} to {}".format(var_name, CONTROL_DESC, TEST2_DESC))
    plt.xlabel("Days since model start")
    plt.ylabel("Pressure (hPa)".format(var_name, var_units))
    plt.colorbar()
    plt.savefig(diff2_filename)
    plt.close()
