"""
Example script for application on ukcp18 data.
"""
import xarray as xr
from eto import ETo_xr
from eto.util_xr import prepare_input_from_ukcp18, get_ukcp18_kwargs
import sys
import glob
import matplotlib.pyplot as plt
from datetime import datetime

varlist = ['hurs', 'psl', 'rls', 'rss', 'sfcWind', 'tasmax', 'tasmin', 'tas']
pd = dict([[v,
       "/Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc".format(v)]
      for v in varlist])
outvarname = 'evpot-fao'
outpath = "/Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc".format(outvarname)

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

et_ukcp18_fao.mean(dim=['projection_x_coordinate', 'projection_y_coordinate']).isel(bnds=0)\
    .groupby('time.year').sum().plot()
plt.show()

# increasing ETo looks to be the result of:
#  decreasing RH_mean (+)
#  increasing R_n (R_s increases more than R_l) (+)
#  increasing T (+)
#  ?? effect of higher P after 2040 or so??
#  slightly decreasing U_z (-)
