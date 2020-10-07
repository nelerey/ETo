"""
Testing the ET0 package Sarah sent.
This computes FAO56 PM well-watered grass PE.
"""
from eto import ETo_xr, datasets
from eto.util_xr import get_ukcp18_kwargs, prepare_input_from_ukcp18
import pandas as pd
import numpy as np

# test on generic dataset.
# this dataset contains:
#   R_sm
#   T_max
#   T_min
#   e_a
# along with date.

et1 = ETo_xr()
ex1_path = datasets.get_path('example_daily')
dftemp = pd.read_csv(ex1_path, parse_dates=True, infer_datetime_format=True, index_col='date')
ds = dftemp.to_xarray().assign_coords({'projection_x_coordinate': 20000,
                                       'projection_y_coordinate': 50000}
                                      ).expand_dims(['projection_x_coordinate',
                                                     'projection_y_coordinate'])

z_msl = 500  # known Elevation of the met station above mean sea level (m) (only needed if P is not in df).
lat = -43.6  # knwon The latitude of the met station (dec deg) (only needed if R_s or R_n are not in df).
lon = 172  # known  The longitude of the met station (dec deg) (only needed if calculating ETo hourly)
TZ_lon = 173  # known The longitude of the center of the time zone (dec deg) (only needed if calculating ETo hourly).
# z_u = 2:  The height of the wind speed measurement (m). Default is 2 m.
freq = 'D'  # known
# et1.param_est(tsdata, freq, z_msl, lat, lon, TZ_lon)
# et1.ts_param.head()
et1.param_est_xr(ds, freq, z_msl, lat, lon, TZ_lon)
et1.ts_param


# test on ukcp18
varlist = ['hurs', 'psl', 'rls', 'rss', 'sfcWind', 'tasmax', 'tasmin', 'tas']
pd = dict([[v,
       "/Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc".format(v)]
      for v in varlist])
ds_ukcp18 = prepare_input_from_ukcp18(pd, concat_dim=False)
kwargs_ukcp18 = get_ukcp18_kwargs()
et_ukcp18 = ETo_xr()
et_ukcp18.param_est_xr(ds_ukcp18, **kwargs_ukcp18)
