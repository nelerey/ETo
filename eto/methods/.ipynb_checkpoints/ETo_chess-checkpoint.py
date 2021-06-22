# -*- coding: utf-8 -*-
"""
calculates potential evaporation from UKCP18 variables based on the method by robinson 2017.
"""
import xarray as xr
import numpy as np
import argparse
import sys
from datetime import datetime


class constants:

    def __init__(self):
        self.t_steam = 373.15  # steam point temperature in kelvin!
        self.p_steam = 101325  # steam point pressure in Pa
        self.a_factors = [13.3185, -1.9760, -0.6445, -0.1299]
        self.epsilon = 0.622  # robinson 2017 eq 3; gill 1982
        self.la = 2.5 * 10 ** 6  # J kg -1; latent heat of evaporation of water
        self.c_p = 1010  # J kg-1 K-1  # specific heat capacity of air
        self.psychrometric = 0.0004  # K-1; psychrometric constant
        # everything is FAO56-grass:
        self.r_s = 70.0  # s m-1  stomatal resistance
        self.alpha = 0.23
        self.emissivity = 0.92
        self.r = 287.05  # J kg-1 K-1, gas consstant of air
        self.t_d = 60*60*24  # seconds day-1: lengthof a day in seconds


def faster_polyval(p, x):
    """
    thank you stackoverflow user Joe Kingston: https://stackoverflow.com/questions/32526388/a-faster-numpy-polynomial
    p=coefficients
    x=array to evalate polynomial
    """
    y = np.zeros(x.shape, dtype=float)
    for i, v in enumerate(p[::-1]):
        y *= x
        y += v
    return y


def richards_poly(t_a, derive=False):
    """
    Computes the polynomial that forms the exponent in eq 5 Robinson et al 2017 (richards approximation of q_s=f(t_a), 1971). if derive=True, it will compute the first derivative of this term to t_a instead
    deze functie geeft het gewenste resultaat.

    Parameters
    ----------
    t_a: np.ndarray of air temperature values in Kelvin

    Return
    ------
    a numpy.ndarray same shape as t_a which is the exponent for eq5
    """

    tterm = (1. - cst.t_steam/t_a)

    if derive:
        print("calculating richards poly derivative")
        coefficients = cst.a_factors * np.arange(1,5)
        res = faster_polyval(coefficients, tterm)
        return res * cst.t_steam / (t_a ** 2)
    else:
        print("calculating richards poly")
        coefficients = np.array([0, *cst.a_factors])
        res = faster_polyval(coefficients, tterm)
        return res


def calc_sat_spechum(e, p_air):
    """
    eq 3 in robinson et al 2017
    """
    q = cst.epsilon * e / (p_air - (1-cst.epsilon) * e)
    return q


def calc_delta_qs(e_s, t_a, q_s, p_air):
    """
    eq6 robinson et al 2017
    """
    return richards_poly(t_a, derive=True) * p_air * q_s / (p_air - (1-cst.epsilon) * e_s)


def calc_delta_es(t_a):
    """
    Derive e_s with respect to t_air.

    Derivation of the implementation:
    d(e_s)/d(t_a) = d(p_steam * exp(richards-polynomial))/d(t_a)
                  = p_steam * exp(richards-polynomial) * d(richards-polynomial)/d(t_a)
    """
    return cst.p_steam * np.exp(richards_poly(t_a)) * richards_poly(t_a, derive=True)


def calc_eto_ceh_xr(t_a, p_air, netlongwavedown, netshortwavedown, u_10, q_a):
    """
    Function to estimate reference ET (ETo) roughly following the method used to produce CHESS-PE ([2]) using the reference crop defined in `FAO 56 paper <http://www.fao.org/docrep/X0490E/X0490E00.htm>`_ [1]_

    I did not calculate the specific humidiy and instead used the vapour pressure deficit as in Monteith et al. 1965, as I don't know why they did not also do that. They refer to some other publication where q_s is used instead of e_s, but there it is also not explained why. This impacts the deficit term in the aerodynamic component and the delta term.

    Parameters
    ----------

    Returns
    -------

    References
    ----------

    .. [1] Allen, R. G., Pereira, L. S., Raes, D., & Smith, M. (1998). Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56. FAO, Rome, 300(9), D05109.

    .. [2] Robinson, E. L., Blyth, E. M., Clark, D. B., Finch, J., & Rudd, A. C. (2017). Trends in atmospheric evaporative demand in Great Britain using high-resolution meteorological data. Hydrology and Earth System Sciences, 21(2), 1189-1224.
    """
    # calculate the derivative of saturated vapour pressure with respect to air temperature. this does not equate eq 6 in Robinson et al. 2017, see function description.
    e_s = cst.p_steam * np.exp(richards_poly(t_a))
    q_s = calc_sat_spechum(e_s, p_air)
    delta_qs = calc_delta_qs(e_s, t_a, q_s, p_air)

    # eq 7: calculate net radiation from the short- and longwave downwelling radiation, with surface temperature being approximated by air temperature and assuming G is zero
    # R_n = SWnet + LWnet = Sd - cst.alpha*Sd + cst.emissivity * Ld - cst.emissivity * cst.stefboltz * T_a**4

    A = netshortwavedown + netlongwavedown  # = R_n - G, but G is assumed to be 0

    # eq 9: air density
    rho_air = p_air / (cst.r*t_a)

    # eq 10
    r_a = 278./u_10  # aerodynamic resistance for a reference crop with height 0.12 m

    # denominator for eq 11/12.
    denominator = delta_qs + cst.psychrometric * (1 + (cst.r_s/r_a))

    # eq 11
    radiative_component = (delta_qs * A) / denominator
    # eq 12
    aerodynamic_component = (cst.c_p*rho_air/r_a) * (q_s - q_a) / denominator
    pe = (radiative_component + aerodynamic_component) * cst.t_d / cst.la
    return pe


def main():
    """

    """
    print("Running main. Start parsing arguments...")

    parser = argparse.ArgumentParser(description="Compute PET using the Robinson et al. 2017 method")
    parser.add_argument('-v', '--variables', action='store',
                        dest='varlist', default='',
                        help=('comma-separated list of variable names to be used in PE calculation. \n'
                              'Expample: huss,psl,rls,rss,sfcWind,tas'))
    parser.add_argument('-i', '--input-path-format', action='store',
                        dest='in_path_format', default='',
                        help=('Generic path to the input files, with {var} where the variable name should go.'
                              'Note that this script assumes this path format to apply to all variables. \n'
                              'Example: /Users/nelereyniers/data/toydata/toydata_3x3/{}_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc'))
    parser.add_argument('-o', '--output-path', action='store',
                        dest='output_path', default='pet.nc',
                        help=('Path for the output file.  '
                              'Example: /Users/nelereyniers/data/toydata/toydata_3x3/pet_rcp85_land-rcm_uk_12km_01_day_19801201-20801130.nc'))

    args = parser.parse_args()
    print("args: \n", args)
    varlist = args.varlist.split(',')
    input_path_format = args.in_path_format
    output_path = args.output_path

    # load data
    pathsdict = dict([[v, input_path_format.format(var=v)] for v in varlist])
    datadict = {}
    for key in pathsdict:
        # if concat_dim:
        #    datadict[key] = xr.open_mfdataset(pathsdict[key], concat_dim=concat_dim,
        #                                      chunks={'projection_x_coordinate': 15, 'projection_y_coordinate': 15})
        #else:
        datadict[key] = xr.open_dataset(pathsdict[key])[key]

    # convert units where needed
    datadict['psl'] = datadict['psl'] * 100  # hPa to Pa
    datadict['tas'] = datadict['tas'] + 273.15  # degC to K
    print("Data loaded and units of psl and tas converted. \ndatadict: \n", datadict)


    pet_values = calc_eto_ceh_xr(
        t_a=datadict['tas'].values,
        p_air=datadict['psl'].values,
        netlongwavedown=datadict['rls'].values,
        netshortwavedown=datadict['rss'].values,
        u_10=datadict['sfcWind'].values,
        q_a=datadict['huss'].values
    )
    pet = xr.DataArray(pet_values, datadict['tas'].coords, datadict['tas'].dims)
    pet.name = "pet"
    print("PE calculated: \n", pet, "\n, start preparing for export...")


    allow_neg=False
    if not allow_neg:
        pet = pet.where(pet > 0., 0)  # pe is preserved where condition is true, "other" is applied where condition is false
        print("removed negatives")

    # prepare for export
    attributes = {
        'creation_script': sys.argv[0] + " (N. Reyniers)",
        'creation_time': datetime.now().strftime("%d-%m-%Y %H:%M"),
        'standard_name': 'water_reference_evaporation_amount',
        'long_name': 'FAO56/Robinson reference evaporation amount',
        'units': 'kg m-2',
        'label_units': 'kg m-2',
        'plot_label': 'ET0 (kg m-2)',
        'description': 'Reference crop potential evaporation (kg m-2) calculated from rls, rss, psl, huss, tas, sfcWind based on Robinson et al. (2017) to match CHESS-PE',
        'units_equivalent': 'mm',
        'PE_calculation_method': 'Penman-Monteith as implemented in Robinson et al. (2017), using the FAO56 reference crop (Allen et al., 1998)',
        'references': '[1] Robinson, E. L., et al. "Trends in atmospheric evaporative demand in Great Britain using high-resolution meteorological data." Hydrology and Earth System Sciences 21.2 (2017): 1189-1224. [2] Allen, R. G., et al. \"Crop evapotranspiration-Guidelines for computing crop water requirements-FAO Irrigation and drainage paper 56.\" Fao, Rome 300.9 (1998): D05109.',
        'input_variable_names': ', '.join(list(pathsdict.keys()))
    }
    input_attributes = dict(
        [['input_{}'.format(vn), pathsdict[vn]] for vn in list(pathsdict.keys())]
    )
    pet = pet.assign_attrs(attributes)
    pet = pet.assign_attrs(input_attributes)
    pet.to_netcdf(output_path, engine='netcdf4', format='NETCDF4',
                  encoding={"pet": {"dtype": "float32"}})
    print("Succesfully exported.")


if __name__ == '__main__':
    cst = constants()
    main()
else:
    print("This file was not run as main FYI")
