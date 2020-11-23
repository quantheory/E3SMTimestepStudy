#!/usr/bin/env python

import os

# Note: uses genutil from UV-CDAT.
import numpy as np
import genutil
import matplotlib
import matplotlib.pyplot as plt

import utils
import core_parameter
from taylor_diagram import TaylorDiagram

# These functions taken from metrics in E3SM diags
def corr(model, obs, axis='xy'):
    corr = -np.infty
    try:
        corr = float(genutil.statistics.correlation(
            model, obs, axis=axis, weights='generate'))
    except Exception as err:
        print(err)

    return corr

def std(variable, axis='xy'):
    std = -np.infty
    try:
        std = float(genutil.statistics.std(
            variable, axis=axis, weights='generate'))
    except Exception as err:
        print(err)

    return std

parameter = core_parameter.CoreParameter()
parameter.test_data_path = '/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon/'
parameter.reference_data_path = '/p/lscratchh/santos36/timestep_monthly_avgs_lat_lon/'
parameter.ref_name = 'timestep_ctrl'
parameter.results_dir = '/g/g14/santos36/Timesteps/'
parameter.seasons = ['ANN']
parameter.variables = ['PRECT', 'PRECL', 'PRECC', 'PSL', 'SWCF', 'LWCF',
                       'TGCLDLWP', 'TGCLDIWP', 'CLDTOT', 'U10', 'OMEGA500', 'TS']

ref_data = utils.dataset.Dataset(parameter, ref=True)

marker = ['o', 'd', '+', 's', '>', '<', 'v', '^', 'x', 'h', 'X', 'H'] 
color = ['k', 'r', 'g', 'y', 'm']

matplotlib.rcParams.update({'font.size': 20})
fig = plt.figure(figsize=(9,8))
refstd = 1.0
taylordiag = TaylorDiagram(refstd, fig=fig, rect=111, label="REF")
ax = taylordiag._ax

case_colors = {
    'timestep_all_10s': 'k',
    'timestep_ctrl_ne16': 'r',
    'timestep_ctrl_y04-y07': 'g',
}

for season in parameter.seasons:
    # Will not actually work for multiple seasons as written due to having only
    # one figure.
    print('Season: {}'.format(season))

    for test_name in ('timestep_all_10s', 'timestep_ctrl_ne16', 'timestep_ctrl_y04-y07'):
        parameter.test_name = test_name
        test_data = utils.dataset.Dataset(parameter, test=True)

        irow = 0

        for var in parameter.variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var

            mv1 = test_data.get_climo_variable(var, season)
            mv2 = ref_data.get_climo_variable(var, season)

            parameter.viewer_descr[var] = mv1.long_name if hasattr(
                mv1, 'long_name') else 'No long_name attr in test data.'

            region = 'global'
            if var == 'TS':
                region = 'land'
            land_frac = test_data.get_climo_variable('LANDFRAC', season)
            ocean_frac = test_data.get_climo_variable('OCNFRAC', season)

            mv1_domain = utils.general.select_region(region, mv1, land_frac, ocean_frac, parameter)
            mv2_domain = utils.general.select_region(region, mv2, land_frac, ocean_frac, parameter)

            metrics = dict()
            metrics['ref_regrid'] = {
                'std': float(std(mv2_domain))
            }
            metrics['test_regrid'] = {
                'std': float(std(mv1_domain))
            }
            metrics['misc'] = {
                'corr': float(corr(mv1_domain, mv2_domain))
            }

            std_norm, correlation = metrics['test_regrid']['std']/metrics['ref_regrid']['std'], metrics['misc']['corr']
            taylordiag.add_sample(std_norm, correlation, marker=marker[irow], c=color[0], ms=10,
                                  label=parameter.viewer_descr[var], markerfacecolor='None',
                                  markeredgecolor=case_colors[test_name], linestyle='None')
            irow += 1

        if test_name == 'timestep_all_10s':
            # Add a legend to the figure.
            fig.legend(taylordiag.samplePoints,
                       ["Reference", "Total Precipitation", "Large-scale Precipitation",
                        "Convective Precipitation",
                        "Sea Level Pressure", "Shortwave Cloud Forcing", "Longwave Cloud Forcing",
                        "Liquid Water Path", "Ice Water Path",
                        "Cloud Fraction",
                        "10m wind speed", "500 mb Vertical Velocity", "Land Surface Temperature"],
                       numpoints=1, loc='center right', bbox_to_anchor=(1.0, .75), prop={'size':10})
#    model_text = 'Test Model: ' + parameter.test_name
#    ax.text(0.6, 1, model_text, ha='left', va='center', transform=ax.transAxes, color=color[0], fontsize=12)
#    ax.text(0.6, 0.95, 'Ref. Model: ' + parameter.ref_name, ha='left', va='center', transform=ax.transAxes, color='k', fontsize=12)

#    plt.title(season + ': Spatial Variability', y=1.08)
#    plt.title('Taylor Diagram - Spatial Variability', y=1.08)
    fig.savefig(os.path.join(parameter.results_dir, season + '_metrics_taylor_diag_timesteps.png'))
