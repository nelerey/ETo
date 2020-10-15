"""
Example script for application on ukcp18 data.
"""
import xarray as xr
import sys
import os
print("os.getcwd: \n", os.getcwd())
sys.path.append(os.getcwd())
print("sys.path: \n", sys.path)
# import glob
from datetime import datetime
from optparse import OptionParser
from eto.util_xr import prepare_input_from_ukcp18, get_ukcp18_kwargs
from eto import ETo_xr

parser = OptionParser()
parser.add_option('-v', '--variables', action='store',
                  type='string', dest='varlist', default='',
                  help=('comma-separated list of variable names to be used in PE calculation. \n'
                        'Expample: hurs,psl,rls,rss,sfcWind,tasmax,tasmin,tas'))
parser.add_option('-i', '--input-path-format', action='store',
                  type='string', dest='in_path_format', default='',
                  help=('Generic path to the input files, with {var} where the variable name should go.'
                        'Note that this script assumes this path format to apply to all variables. \n'
                        'Example: /Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc'))
parser.add_option('-o', '--output-path-format', action='store',
                  type='string', dest='out_path_format', default='{}.nc',
                  help=('Path for the output file, where {} will be substituted by the variable name \n'
                       'supplied under the -n (--name) option.'
                        'Example: /Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc'))
parser.add_option('-n', '--name', action='store',
                  type='string', dest='outvarname', default='evpot-result',
                  help='Name to store the output dataset under. \nExample: evpot-fao')

(options, args) = parser.parse_args()
varlist = options.varlist.split(',')
input_path_format = options.in_path_format
output_path_format = options.out_path_format
outvarname = options.outvarname

# varlist = ['hurs', 'psl', 'rls', 'rss', 'sfcWind', 'tasmax', 'tasmin', 'tas']
# pd = dict([[v,
#       "/Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc".format(v)]
#      for v in varlist])
# outvarname = 'evpot-fao'
# outpath = "/Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc".format(outvarname)

# ---------------------------------------------------------------------
pd = dict([[v, input_path_format.format(var=v)] for v in varlist])
outpath = output_path_format.format(outvarname)

ds_ukcp18 = prepare_input_from_ukcp18(pd, concat_dim=False)
kwargs_ukcp18 = get_ukcp18_kwargs()
et_ukcp18 = ETo_xr()
et_ukcp18.param_est_xr(ds_ukcp18, **kwargs_ukcp18)
et_ukcp18_fao = et_ukcp18.eto_fao_xr(remove_extreme_values=False, round_decimals=2)

attributes = {
    'creation_script': sys.argv[0] + " (github user mullenkamp and N. Reyniers)",
    'creation_time': datetime.now().strftime("%d-%m-%Y %H:%M"),
    'standard_name': 'water_potential_evaporation_amount',
    'long_name': 'FAO reference crop evaporation amount',
    'units': 'kg m-2',
    'label_units': 'kg m-2',
    'plot_label': 'FAO56-Penman-Monteith ETo (kg m-2)',
    'description': 'FAO56-Penman-Monteith reference crop potential evaporation (kg m-2)',
    'units_equivalent': 'mm',
    'PE_calculation_method': 'FAO56 Penman-Monteith (Allen et al., 1998)',
    'references': 'Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines \
     for computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300(9), D05109. \
     Procedure implemented from http://www.fao.org/docrep/X0490E/X0490E00.htm',
    'input_variable_names': ', '.join(list(pd.keys()))
}
input_attributes = dict(
    [['input_{}'.format(vn), pd[vn]] for vn in list(pd.keys())]
)
et_ukcp18_fao.assign_attrs(attributes)
et_ukcp18_fao.assign_attrs(input_attributes)
et_ukcp18_fao.name = outvarname
et_ukcp18_fao.to_netcdf(outpath, engine='netcdf4')
