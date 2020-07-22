#!/usr/bin/env python
"""Some remote data access routines to access the NCEO 
data storage on JASMIN."""
import datetime as dt

import gdal
import pandas as pd
import xarray as xr

JASMIN_URL = "http://gws-access.ceda.ac.uk/public/odanceo/"

ERA5_VARIABLES = [
    "hum",
    "precip",
    "ssrd",
    "t2m_max",
    "t2m_mean",
    "t2m_min",
    "wspd",
]
MODIS_VARIABLES = ["Fpar_500m", "Lai_500m", "FparLai_QC"]


def get_era5_ds(variable, remote_url=JASMIN_URL):
    """Return an xarray dataset with a ERA5 variable.
    Function takes care of getting the timestamps
    sorted, and selects all available time steps.

    Parameters
    ----------
    variable: str
        The variable of interest. Must be one of `hum`,
        `precip`, `ssrd`, `t2m_max`, `t2m_mean, `t2m_min
        or `wspd`.
    remote_url: str
        The parent URL. By default, the main NCEO ODA public
        HTML folder
    
    Returns
    -------
    An all singing, all dancing xarray dataset
    """
    assert (
        variable in ERA5_VARIABLES
    ), f"{variable} not one of {ERA5_VARIABLES}"
    today = dt.datetime.now()
    arrays = []
    for year in range(2002, today.year + 1):
        url = f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_{year}.tif"
        retval = gdal.Info(url, allMetadata=True, format="json")
        dates = [
            pd.to_datetime(d["metadata"][""]["Date"]) for d in retval["bands"]
        ]
        ds = xr.open_rasterio(url, chunks={})
        ds = ds.rename({"band": "time"})
        ds = ds.assign_coords({"time": dates})
        arrays.append(ds)

    ds = xr.concat(arrays, dim="time")
    return ds


def get_modis_ds(remote_url=JASMIN_URL, product="Fpar_500m"):
    """Return an xarray dataset with a MODIS variable.
    Function takes care of getting the timestamps
    sorted, and selects all available time steps.

    Parameters
    ----------
    variable: str
        The variable of interest. Must be one of "Fpar_500m",
        `Lai_500m` or `FparLai_QC`
    remote_url: str
        The parent URL. By default, the main NCEO ODA public
        HTML folder
    
    Returns
    -------
    An all singing, all dancing xarray dataset
    TODO: We should really apply the QC layer
    """
    assert (
        product in MODIS_VARIABLES
    ), f"{product} is not one of {MODIS_VARIABLES}"
    today = dt.datetime.now()
    arrays = []
    for year in range(2002, today.year + 1):
        url = f"/vsicurl/{remote_url}/MCD15/{product}_{year}.tif"
        retval = gdal.Info(url, allMetadata=True, format="json")
        dates = [
            pd.to_datetime(
                f"{year}/{d['metadata']['']['DoY']}", format="%Y/%j"
            )
            for d in retval["bands"]
        ]
        ds = xr.open_rasterio(url, chunks={})
        ds = ds.rename({"band": "time"})
        ds = ds.assign_coords({"time": dates})
        arrays.append(ds)

    ds = xr.concat(arrays, dim="time")
    return ds


def calculate_climatology(
    ds, first_year=2002, period="time.month", last_year=None
):
    """A simple function to calculate the mean and standard deviation for
    every `period` (e.g. month or whatever) on a pixel by pixel basis.
    The function allows you to specify a `first_year` and `last_year` for the
    calculation. 

    Parameters
    -----------
    ds: xarray DataFrame
        An xarray dataframe indexed by time. Needs a `time` dimension.
    first_year: int
        The first year of the climatology.
    last_year: int
        The last year of the climatology.
    period: str
        A pandas/xarray time frequency string. Default is monthly.
    """
    if last_year is None:
        last_year = dt.datetime.now().year

    clim_mean = (
        ds.sel(time=slice(f"{first_year}-01-01", f"{last_year}-01-01"))
        .groupby(period)
        .mean("time")
    )
    clim_std = (
        ds.sel(time=slice(f"{first_year}-01-01", f"{last_year}-01-01"))
        .groupby(period)
        .std("time")
    )
    return clim_mean, clim_std
