# -*- coding: utf-8 -*-
"""
Function to estimate reference ET (ETo) from the FAO 56 paper using a minimum of T_min and T_max for daily estimates and T_mean and RH_mean for hourly, but utilizing the maximum number of available met parameters.

xr by nele reyniers: TODO
"""
import numpy as np
import pandas as pd
import os
import xarray as xr
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

def test_eto_fao_daily():
    et1 = ETo_xr(tsdata, 'D', z_msl, lat, lon, TZ_lon, remove_extreme_values=False, round_decimals=2)
    eto1 = et1.eto_fao_xr().sum()
    res1 = tsresults['ETo_FAO_mm'].sum()

    assert eto1 == res1


et1 = ETo_xr(tsdata, 'D', z_msl, lat, lon, TZ_lon, remove_extreme_values=False, round_decimals=2)

def test_eto_har_daily():
    eto2 = et1.eto_hargreaves_xr().sum()
    res1 = tsresults['ETo_Har_mm'].sum()

    assert eto2 == res1


def test_eto_fao_hourly():
    tsdata2 = et1.ts_param_xr[['R_s', 'T_mean', 'e_a']]
    tsdata3 = et1.tsreg(tsdata2, 'H', 'time')
    et2 = ETo_xr(tsdata3, 'H', z_msl, lat, lon, TZ_lon, remove_extreme_values=False, round_decimals=2)
    eto3 = et2.eto_fao_xr().sum()

    res1 = tsresults['ETo_FAO_mm'].sum()

    assert eto3 > res1
