"""
Utility functions for the UKCP18-based thing
"""
import xarray as xr

def prepare_input_from_ukcp18(pathsdict, concat_dim='time'):
    """
    Function to prepare a xr.Dataset to use as direct input for the eto_xr functions and par_est_xr.
    This is developed for use on data from one

    Parameters
    ----------
    pathsdict: (dict)
     dictionary with the ukcp18 short variable names as keys and the paths or lists of paths as items.
    concat_dim: (str)
     to be supplied to xr.open_mfdataset. use concat_dim=False in case xr.open_dataset should be used. TODO

    Returns
    -------
    xr.Dataset to use as direct input for the eto_xr functions and par_est_xr
    """
    # open the datasets
    datadict = {}
    for key in pathsdict:
        if concat_dim:
            datadict[key] = xr.open_mfdataset(pathsdict[key], concat_dim=concat_dim)
        else:
            datadict[key] = xr.open_dataset(pathsdict[key])
    # put ds containing different variables into one dataset
    ds_ukcp18 = xr.merge(list(datadict.values()))
    ds_ukcp18 = ds_ukcp18.rename_vars({'hurs': 'RH_mean',
                                       'tas': 'T_mean',
                                       'tasmin': 'T_min',
                                       'tasmax': 'T_max',
                                       'sfcWind': 'U_z',
                                       'psl': 'P',
                                       'rls': 'R_nl',
                                       'rss': 'R_ns'})
    # unit conversions
    ds_ukcp18['P'] = ds_ukcp18['P']/10  # hPa to kPa
    ds_ukcp18['R_nl'] = ds_ukcp18['R_nl'] * -0.0864  # Wm^-2 to MJm^-2, daily resolution, upwelling
    ds_ukcp18['R_ns'] = ds_ukcp18['R_ns'] * 0.0864  # Wm^-2 to MJm^-2, daily resolution

    return ds_ukcp18


def get_ukcp18_kwargs():
    """
    From UKCP18 metadata.

    Returns
    -------
    dict with UKCP18-appropriate values of keyword arguments for par_est_xr and eto_xr method functions.
    """
    ukcp18rcm_kwargs = {
        'freq': 'D',
        'z_u': 10,
        'dt_index_name': 'time',
    }
    return ukcp18rcm_kwargs
