#!/usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import basemap
import netCDF4 as nc4

cmap = plt.get_cmap('coolwarm')
bmap = basemap.Basemap(lon_0=180.)

suffix = '_avg_d01-10'

REF_FILE_NAMES = []
TEST_FILE_NAMES = []

for day in ["01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]:
    for time in ["00000", "03600", "07200", "10800", "14400", "18000", "21600",
                 "25200", "28800", "32400", "36000", "39600", "43200", "46800",
                 "50400", "54000", "57600", "61200", "64800", "68400", "72000",
                 "75600", "79200", "82800"]:
        if day == "01" and time == "00000":
            continue
        REF_FILE_NAMES.append('/p/lscratchh/santos36/interpolated/timestep_ctrl_atm_01-{}-{}.nc'.format(day, time))
        TEST_FILE_NAMES.append('/p/lscratchh/santos36/interpolated/timestep_CLUBB_MG2_10s_atm_01-{}-{}.nc'.format(day, time))

REF_FILE_NAMES.append('/p/lscratchh/santos36/interpolated/timestep_ctrl_atm_01-06-00000.nc')
TEST_FILE_NAMES.append('/p/lscratchh/santos36/interpolated/timestep_CLUBB_MG2_10s_atm_01-06-00000.nc')

rfiles = [nc4.Dataset(ref_name, 'r') for ref_name in REF_FILE_NAMES]
tfiles = [nc4.Dataset(test_name, 'r') for test_name in TEST_FILE_NAMES]

num_files = len(REF_FILE_NAMES)
assert num_files == len(TEST_FILE_NAMES)
print("number of files is: ", num_files)

lat = rfiles[0]['lat']
lon = rfiles[0]['lon']
lev = rfiles[0]['lev']
ilev = rfiles[0]['ilev']
gw = rfiles[0]['gw']

print("Sum of gw:", gw[:].sum())
gw = gw[:]/gw[:].sum()

ref_swcf = sum([rfile['SWCF'][0,:,:] for rfile in rfiles]) / num_files
test_swcf = sum([tfile['SWCF'][0,:,:] for tfile in tfiles]) / num_files

diff_swcf = test_swcf - ref_swcf
diff_swcf_weighted = diff_swcf.copy()
diff_swcf2_weighted = diff_swcf**2
for i in range(len(lat)):
    diff_swcf_weighted[i,:] *= gw[i]
    diff_swcf2_weighted[i,:] *= gw[i]

swcf_mean = diff_swcf_weighted.sum()/len(lon)
swcf_rmse = np.sqrt(diff_swcf2_weighted.sum()/len(lon))

print("SWCF mean diff", swcf_mean)
print("SWCF RMSE: ", swcf_rmse)

plt.pcolormesh(lon[:], lat[:], ref_swcf)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([-850.,0.])
plt.savefig('SWCF_ctrl{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], test_swcf)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([-850.,0.])
plt.savefig('SWCF_CLUBB_MG2{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], diff_swcf, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([-3., 3.])
plt.savefig('SWCF_diff_CLUBB_MG2{}.png'.format(suffix))
plt.close()

ref_precc = sum([rfile['PRECC'][0,:,:] for rfile in rfiles]) / num_files
test_precc = sum([tfile['PRECC'][0,:,:] for tfile in tfiles]) / num_files

ref_prect = ref_precc + sum([rfile['PRECL'][0,:,:] for rfile in rfiles]) / num_files
test_prect = ref_prect + sum([tfile['PRECL'][0,:,:] for tfile in tfiles]) / num_files

ref_precc *= 1000. * 86400.
test_precc *= 1000. * 86400.
ref_prect *= 1000. * 86400.
test_prect *= 1000. * 86400.

diff_precc = test_precc - ref_precc
diff_precc_weighted = diff_precc.copy()
diff_precc2_weighted = diff_precc**2
for i in range(len(lat)):
    diff_precc_weighted[i,:] *= gw[i]
    diff_precc2_weighted[i,:] *= gw[i]

precc_mean = diff_precc_weighted.sum()/len(lon)
precc_rmse = np.sqrt(diff_precc2_weighted.sum()/len(lon))

print("PRECC mean diff", precc_mean)
print("PRECC RMSE: ", precc_rmse)

diff_prect = test_prect - ref_prect
diff_prect_weighted = diff_prect.copy()
diff_prect2_weighted = diff_prect**2
for i in range(len(lat)):
    diff_prect_weighted[i,:] *= gw[i]
    diff_prect2_weighted[i,:] *= gw[i]

prect_mean = diff_prect_weighted.sum()/len(lon)
prect_rmse = np.sqrt(diff_prect2_weighted.sum()/len(lon))

print("PRECT mean diff", prect_mean)
print("PRECT RMSE: ", prect_rmse)

plt.pcolormesh(lon[:], lat[:], ref_precc)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([0.,65.])
plt.savefig('PRECC_ctrl{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], test_precc)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([0.,65.])
plt.savefig('PRECC_CLUBB_MG2{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], diff_precc, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([-3., 3.])
plt.savefig('PRECC_diff_CLUBB_MG2{}.png'.format(suffix))
plt.close()



plt.pcolormesh(lon[:], lat[:], ref_prect)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([0.,35.])
plt.savefig('PRECT_ctrl{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], test_prect)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([0.,35.])
plt.savefig('PRECT_CLUBB_MG2{}.png'.format(suffix))
plt.close()

plt.pcolormesh(lon[:], lat[:], diff_prect, cmap=cmap)
bmap.drawcoastlines()
ax = plt.gca()
plt.axis('tight')
ax.set_xticks([0., 90., 180., 270., 360.])
ax.set_xticklabels(['0', '90E', '180', '90W', '0'])
ax.set_yticks([60., 30., 0., -30., -60.])
ax.set_yticklabels(['60N', '30N', '0', '30S', '60S'])
plt.colorbar()
#plt.clim([-3., 3.])
plt.savefig('PRECT_diff_CLUBB_MG2{}.png'.format(suffix))
plt.close()
