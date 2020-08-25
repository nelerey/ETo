
from eto import ETo, datasets
import numpy as np
import pandas as pd
import xarray as xr


self = ETo()
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


#def param_est(self, df, freq='D', z_msl=None, lat=None, lon=None, TZ_lon=None,
#              z_u=2, K_rs=0.16, a_s=0.25, b_s=0.5, alb=0.23, dt_index_name='date):
z_u = 2
K_rs = 0.16
a_s = 0.25
b_s = 0.5
alb = 0.23
dt_index_name = 'date'

"""
Function to estimate the parameters necessary to calculate reference ET (ETo) from the `FAO 56 paper
<http://www.fao.org/docrep/X0490E/X0490E00.htm>`_ [1]_ using a minimum of T_min and T_max for daily
estimates and T_mean and RH_mean for hourly, but optionally utilising the maximum number of available
met parameters. The function prioritizes the estimation of specific parameters based on the available input data.

Parameters
----------
df : DataFrame
    Input Metereological data (see Notes section).
z_msl : float, int, or None
    Elevation of the met station above mean sea level (m) (only needed if P is not in df).
lat : float, int, or None
    The latitude of the met station (dec deg) (only needed if R_s or R_n are not in df).
lon : float, int, or None
    The longitude of the met station (dec deg) (only needed if calculating ETo hourly)
TZ_lon : float, int, or None
    The longitude of the center of the time zone (dec deg) (only needed if calculating ETo hourly).
z_u : float or int
    The height of the wind speed measurement (m). Default is 2 m.
freq : str
    The Pandas time frequency string of the input and output. The minimum frequency is hours (H) and the maximum is month (M).
K_rs : float
    Rs calc coefficient (0.16 for inland stations, 0.19 for coastal stations)
a_s : float
    Rs calc coefficient
b_s : float
    Rs calc coefficient
alb : float
    Albedo. Should be 0.23 for the reference crop.

Returns
-------
DataFrame

Notes
--------
The input data must be a DataFrame with specific column names according to the met parameter. The column names
should be a minimum of T_min and T_max for daily estimates and T_mean and RH_mean for hourly,
but can contain any/all of the following:

R_n
    Net radiation (MJ/m2)
R_s
    Incoming shortwave radiation (MJ/m2)
G
    Net soil heat flux (MJ/m2)
T_min
    Minimum Temperature (deg C)
T_max
    Maximum Temperature (deg C)
T_mean
    Mean Temperature (deg C)
T_dew
    Dew point temperature (deg C)
RH_min
    Minimum relative humidity
RH_max
    Maximum relative humidity
RH_mean
    Mean relative humidity
n_sun
    Number of sunshine hours per day
U_z
    Wind speed at height z (m/s)
P
    Atmospheric pressure (kPa)
e_a
    Actual Vapour pressure derrived from RH


Parameter estimation values refer to the quality level of the input parameters into the ETo equations.
Where a 0 (or nothing) refers to no necessary parameter estimation (all measurement data was available),
while a 1 refers to parameters that have the best input estimations and up to a value of 3 is the worst.
Starting from the right, the first value refers to U_z, the second value refers to G, the third value
refers to R_n, the fourth value refers to R_s, the fifth value refers to e_a, the sixth value refers
 to T_mean, the seventh value refers to P.

References
----------

.. [1] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines
 for computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300(9), D05109.
"""

met_names = np.array(['R_n', 'R_s', 'G', 'T_min', 'T_max', 'T_mean', 'T_dew',
                      'RH_min', 'RH_max', 'RH_mean', 'n_sun', 'U_z', 'P', 'e_a'])
self.freq = freq
varnames = list(ds.keys())
####################################
##### Set up the DataFrame and estimated values series
new_cols = met_names[~np.in1d(met_names, varnames)]
new_ds = xr.Dataset(dict([[nc, xr.DataArray(np.full(tuple(ds.dims.values()), np.nan),
                                            coords=ds.coords,
                                            dims=ds.dims,
                                            name=nc)]
                          for nc in new_cols]))
self.ts_param = xr.merge([ds, new_ds]).copy()
xr.merge([new_ds, ds])

self.est_val = xr.DataArray(np.full(tuple(ds.dims.values()), 0),
                            coords=ds.coords,
                            dims=ds.dims,
                            name='est_val')

####################################
###### Check to make sure minimum requirements are met
if 'H' in freq:
    T_mean_bool = self.ts_param['T_mean'].isnull().any()
    RH_mean_bool = self.ts_param['RH_mean'].isnull().any()
    e_a_bool = self.ts_param['e_a'].isnull().any()
    if T_mean_bool | (RH_mean_bool & e_a_bool):
        raise ValueError('Minimum data input (T_mean and RH_mean or e_a) was not met. Check your data.')
else:
    T_min_bool = self.ts_param['T_min'].isnull().any()
    T_max_bool = self.ts_param['T_max'].isnull().any()
    if T_min_bool | T_max_bool:
        raise ValueError('Minimum data input (T_min and T_max) was not met. Check your data.')


####################################
###### Calculations

######
### Time index
# if type(df.index) is not pd.DatetimeIndex:
#     raise ValueError('DataFrame must have a datetime index!')

if type(ds.indexes[dt_index_name]) is not pd.DatetimeIndex:
    raise ValueError('DataFrame must have a datetime index!')

# Create the Day of the year vector
Day = ds.indexes[dt_index_name].dayofyear

######
## Atmospheric components

# Air Pressure
# nnn: missing data is filled in by a function of elevation above sea level in ts_param.
# nnn: the value where P is missing in est_val increases by 1 million
# self.est_val.loc[self.ts_param['P'].isnull()] += 1000000
# self.ts_param.loc[self.ts_param.isnull(), 'P'] = 101.3*((293 - 0.0065*z_msl)/293)**5.26

self.est_val = self.est_val + self.ts_param['P'].isnull()*1000000
self.ts_param['P'].where(self.ts_param['P'].notnull(), drop=False,
                         other=101.3*((293 - 0.0065*z_msl)/293)**5.26)

# Psychrometric constant
self.ts_param['gamma'] = (0.665*10**-3)*self.ts_param['P']

# ######
## Temperature and humidity components
self.est_val = self.est_val + self.ts_param['T_mean'].isnull()*1000000

self.ts_param['T_mean'] = xr.where(self.ts_param['T_mean'].isnull(),  # condition
                                   (self.ts_param['T_max'] + self.ts_param['T_min']) / 2,  # value if True
                                   self.ts_param['T_mean'])   # value if False

# ## Vapor pressures
if 'H' in freq:
    # self.ts_param.loc[self.ts_param['e_a'].isnull(), 'e_a'] = \
    #     self.ts_param.loc[self.ts_param['e_a'].isnull(), 'e_mean'] * \
    #     self.ts_param.loc[self.ts_param['e_a'].isnull(), 'RH_mean']/100
    self.ts_param['e_mean'] = 0.6108 * np.exp(17.27 * self.ts_param['T_mean'] /
                                              (self.ts_param['T_mean'] + 237.3))
    self.ts_param['e_a'] = xr.where(self.ts_param['e_a'].isnull(),  # condition
                                   self.ts_param['e_mean'] * self.ts_param['RH_mean'] / 100,  # value if True
                                   self.ts_param['e_a'])   # value if False

else:
    self.ts_param['e_max'] = 0.6108 * np.exp(17.27 * self.ts_param['T_max'] /
                                             (self.ts_param['T_max'] + 237.3))
    self.ts_param['e_min'] = 0.6108 * np.exp(17.27 * self.ts_param['T_min'] /
                                           (self.ts_param['T_min'] + 237.3))
    self.ts_param['e_s'] = (self.ts_param['e_max'] + self.ts_param['e_min']) / 2

    # e_a where dewpoint temperature is known but e_a is not known
    self.ts_param['e_a'] = xr.where(self.ts_param['e_a'].isnull(),
                                    0.6108 * np.exp(17.27 * self.ts_param['T_dew'] /
                                                    (self.ts_param['T_dew'] + 237.3)),
                                    self.ts_param['e_a'])

    # e_a where min and max temperatures and humidities are known, but T_dew and e_a are not known
    self.est_val = self.est_val + self.ts_param['T_dew'].isnull()*10000
    self.ts_param['e_a'] = xr.where(self.ts_param['e_a'].isnull(),
                                    (self.ts_param['e_min'] * self.ts_param['RH_max']/100 +
                                     self.ts_param['e_max'] * self.ts_param['RH_min']/100) / 2,
                                    self.ts_param['e_a'])  # TODO: klopt het dat rhmaxmin en emaxmin gekruist zijn???

    # self.ts_param['e_a'] if only mean humidity is known
    self.est_val = self.est_val + self.ts_param['e_a'].isnull()*10000

    self.ts_param['e_a'] = xr.where(self.ts_param['e_a'].isnull(),
                                    self.ts_param['RH_mean']/100 *
                                    (self.ts_param['e_max'] + self.ts_param['e_min'])/2,
                                    self.ts_param['e_a'])

    # e_a if humidity is not known: use T_min to approximate T_dew
    self.est_val = self.est_val + self.ts_param['e_a'].isnull()*10000  # voeg waarde bij locaties die nu nog 0 zijn

    self.ts_param['e_a'] = xr.where(self.ts_param['e_a'].isnull(),
                                    0.6108 * np.exp(17.27 * self.ts_param['T_min'] /
                                                    (self.ts_param['T_min'] + 237.3)),
                                    self.ts_param['e_a'])

# Delta
self.ts_param['delta'] = 4098 * (0.6108 * np.exp(17.27 * self.ts_param['T_mean'] /
                                                 (self.ts_param['T_mean'] + 237.3))) \
                         / ((self.ts_param['T_mean'] + 237.3) ** 2)


######
## Radiation components

# R_a
phi = lat * np.pi/180  # TODO: waar wordt lat gedefined? want dit moet komen van de coordinaten...
delta = 0.409 * np.sin(2 * np.pi * Day/365-1.39)
d_r = 1 + 0.033 * np.cos(2 * np.pi * Day/365)  # TODO: make 360Day calender compatible
w_s = np.arccos(-np.tan(phi) * np.tan(delta))


if 'H' in freq:
    hour_vec = ds.indexes[dt_index_name].hour
    b = (2 * np.pi * (Day - 81)) / 364
    S_c = 0.1645 * np.sin(2 * b) - 0.1255 * np.cos(b) - 0.025 * np.sin(b)
    w = np.pi/12 * (((hour_vec + 0.5) + 0.6666667 * (TZ_lon - lon) + S_c) - 12)
    w_1 = w - (np.pi * 1)/24  # Need to update one day for different hourly periods
    w_2 = w + (np.pi * 1)/24  # Need to update one day for different hourly periods
    self.ts_param['R_a'] = 12 * 60 / np.pi * 0.082 * d_r * \
                           ((w_2 - w_1) * np.sin(phi) * np.sin(delta) +
                            np.cos(phi) * np.cos(delta) * (np.sin(w_2) - np.sin(w_1)))
else:
    self.ts_param['R_a'] = 24 * 60 / np.pi * 0.082 * d_r * \
                           (w_s * np.sin(phi) * np.sin(delta) +
                            np.cos(phi) * np.cos(delta) * np.sin(w_s))


# Daylight hours
N = 24 * w_s / np.pi  # TODO: make 3D for R_s calculation? what is this

# R_s if n_sun is known
self.est_val = self.est_val + self.ts_param['R_s'].isnull()*1000
self.ts_param['R_s'] = xr.where(self.ts_param['R_s'].isnull,
                                (a_s + b_s * np.squeeze(self.ts_param['n_sun'].values) / N) * self.ts_param['R_a'],
                                self.ts_param['R_s'])  # np.squeeze: necessary when some dimensions have length of 1

# R_s if n_sun is not known
self.est_val = self.est_val + self.ts_param['R_s'].isnull()*1000
self.ts_param['R_s'] = xr.where(self.ts_param['R_s'].isnull(),
                                K_rs * ((self.ts_param['T_max'] - self.ts_param['T_min'])**0.5) * self.ts_param['R_a'],
                                self.ts_param['R_s'])

# R_so
R_so = (0.75 + 2 * 10**(-5) * z_msl) * self.ts_param['R_a']

# R_ns from R_s
R_ns = (1 - alb) * self.ts_param['R_s']

# # R_nl
if 'H' in freq:
    R_nl = (2.043 * 10**(-10)) * ((self.ts_param['T_mean'] + 273.16)**4) * \
           (0.34 - 0.14 * (self.ts_param['e_a'])**0.5) * \
           ((1.35 * self.ts_param['R_s'] / R_so) - 0.35)
else:
    R_nl = (4.903 * 10**(-9)) * (((self.ts_param['T_max'] + 273.16) ** 4 +
                                  (self.ts_param['T_min'] + 273.16) ** 4) / 2) * \
           (0.34 - 0.14 * (self.ts_param['e_a']) ** 0.5) * \
           ((1.35 * self.ts_param['R_s']/R_so) - 0.35)

# R_n
self.est_val = self.est_val + self.ts_param['R_n'].isnull()*100
self.ts_param['R_n'] = xr.where(self.ts_param['R_n'].isnull(),
                                R_ns - R_nl, self.ts_param['R_n'])

# G
self.est_val = self.est_val + self.ts_param['G'].isnull()*10
self.ts_param['G'] = xr.where(self.ts_param['G'].isnull(),
                              0, self.ts_param['G'])

######
## Wind component

self.ts_param['U_2'] = self.ts_param['U_z'] * 4.87 / np.log(67.8 * z_u - 5.42)

# or use 2 if wind speed is not known
self.est_val = self.est_val + self.ts_param['U_z'].isnull()*1
self.ts_param['U_2'] = xr.where(self.ts_param['U_z'].isnull(),
                                2, self.ts_param['U_z'])
