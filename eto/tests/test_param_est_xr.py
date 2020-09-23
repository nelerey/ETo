"""
Testing the ET0 package Sarah sent.
This computes FAO56 PM well-watered grass PE.
"""
from eto import ETo, datasets
import pandas as pd
import numpy as np

# test on generic dataset.
# this dataset contains:
#   R_sm
#   T_max
#   T_min
#   e_a
# along with date.

et1 = ETo()
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
