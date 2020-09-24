# -*- coding: utf-8 -*-
"""
Function to estimate reference ET (ETo) from the FAO 56 paper using a minimum of T_min and T_max for daily estimates and T_mean and RH_mean for hourly, but utilizing the maximum number of available met parameters.

xr by nele reyniers: TODO
"""
import numpy as np
import pandas as pd
import os
from eto import ETo_xr

###############################
### Parameters

_module_path = os.path.dirname(__file__)
example_csv = 'example_daily.csv'
results_csv = 'example_daily_results.csv'

example1 = os.path.join(_module_path, example_csv)
results1 = os.path.join(_module_path, results_csv)

z_msl = 500
lat = -43.6
lon = 172
TZ_lon = 173

###############################
### Tests

tsdata = pd.read_csv(example1, parse_dates=True, infer_datetime_format=True, index_col='date')
tsresults = pd.read_csv(results1, parse_dates=True, infer_datetime_format=True, index_col='date')
ds = tsdata.to_xarray().assign_coords({'projection_x_coordinate': 20000,
                                       'projection_y_coordinate': 50000}
                                      ).expand_dims(['projection_x_coordinate',
                                                     'projection_y_coordinate'])

self = ETo_xr(ds, 'D', z_msl, lat, lon, TZ_lon)
# eto1 = et1.eto_fao().sum()
# res1 = tsresults['ETo_FAO_mm'].sum()
print("self:\n", self)
remove_extreme_values = False
round_decimals = 2
max_ETo = 15
min_ETo = 0
interp = False
maxgap = 15
# def eto_fao(self, max_ETo=15, min_ETo=0, interp=False, maxgap=15, remove_extreme_values=True, round_decimals=2):
"""
Function to estimate reference ET (ETo) from the `FAO 56 paper <http://www.fao.org/docrep/X0490E/X0490E00.htm>`_ [1]_ 
using a minimum of T_min and T_max for daily estimates and T_mean and RH_mean for hourly, but optionally utilising 
the maximum number of available met parameters. 
The function prioritizes the estimation of specific parameters based on the available input data.

Parameters
----------
max_ETo : float or int
    The max realistic value of ETo (mm).
min_ETo : float or int
    The min realistic value of ETo (mm).
interp : False or str
    Should missing values be filled by interpolation? Either False if no interpolation should be performed, or a string of the interpolation method. See Pandas interpolate function for methods. Recommended interpolators are 'linear' or 'pchip'.
maxgap : int
    The maximum missing value gap for the interpolation.

Returns
-------
DataFrame or Series
    If fill=False, then the function will return a Series of estimated ETo in mm. If fill is a str, then the function will return a DataFrame with an additional column for the filled ETo value in mm.

References
----------

.. [1] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for
computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300(9), D05109.
"""

######
## ETo equation
if 'H' in self.freq:
    ETo_FAO = (0.408 * self.ts_param['delta'] * (self.ts_param['R_n'] - self.ts_param['G']) +
               self.ts_param['gamma'] * 37/(self.ts_param['T_mean'] + 273) * self.ts_param['U_2'] *
               (self.ts_param['e_mean'] - self.ts_param['e_a'])) / \
              (self.ts_param['delta'] + self.ts_param['gamma']*(1 + 0.34*self.ts_param['U_2']))
else:
    ETo_FAO = (0.408*self.ts_param['delta']*(self.ts_param['R_n'] - self.ts_param['G']) +
               self.ts_param['gamma'] * 900/(self.ts_param['T_mean'] + 273)*self.ts_param['U_2'] *
               (self.ts_param['e_s'] - self.ts_param['e_a'])) / \
              (self.ts_param['delta'] + self.ts_param['gamma'] * (1 + 0.34 * self.ts_param['U_2']))

ETo_FAO.name = 'ETo_FAO_mm'

if remove_extreme_values:
    print('TODO, extreme values were not removed')  # TODO
    ## Remove extreme values
    # ETo_FAO[ETo_FAO > max_ETo] = np.nan  # TODO, no 3D bool indexing
    # ETo_FAO[ETo_FAO < min_ETo] = np.nan  # TODO, no 3D bool indexing

## ETo equation with filled holes using interpolation (use with caution)
# if isinstance(interp, str):
#     ETo_FAO_fill = self.tsreg(ETo_FAO, self.freq, interp, maxgap)
#     ETo_FAO_fill.name = 'ETo_FAO_interp_mm'
#     ETo = pd.concat([ETo_FAO, ETo_FAO_fill], axis=1).round(2)
# else:
#     ETo = ETo_FAO.round(2)
ETo = ETo_FAO.round(round_decimals)
# return ETo
