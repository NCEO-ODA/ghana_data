#!/usr/bin/env python
"""Some remote data access routines to access the NCEO 
data storage on JASMIN."""
import datetime as dt
import subprocess
from concurrent.futures import ThreadPoolExecutor

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import gdal
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import xarray as xr
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from tqdm.auto import tqdm

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


def get_epsg_code(ds_name):
    """Need to get the EPSG code to pass to cartopy. Not everyone uses
    Lat/Long grids, don't you know...
    """
    # Easiest way seems to be able to pipe out to `gdalsrsinfo` and hope
    # it works
    p = subprocess.Popen(
        ["gdalsrsinfo", "-e", "-o", "epsg", ds_name], stdout=subprocess.PIPE
    )

    std_oot, _ = p.communicate()
    # This might fail?
    epsg_code = std_oot.decode("utf-8").split(":")[1].strip()

    return epsg_code


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
    sample_url = f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_2004.tif"
    epsg = get_epsg_code(sample_url)

    def do_one_year(year):
        url = f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_{year}.tif"
        retval = gdal.Info(url, allMetadata=True, format="json")
        dates = [
            pd.to_datetime(d["metadata"][""]["Date"]) for d in retval["bands"]
        ]
        ds = xr.open_rasterio(url, chunks={"x": 256, "y": 256})
        ds = ds.rename({"band": "time"})
        ds = ds.assign_coords({"time": dates})
        return ds

    with ThreadPoolExecutor(max_workers=8) as executor:
        years = [y for y in range(2002, today.year + 1)]
        arrays = list(
            tqdm(executor.map(do_one_year, years), total=len(years))
        )

    ds = xr.concat(arrays, dim="time")
    ds.attrs["epsg"] = epsg
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
    sample_url = f"/vsicurl/{remote_url}/MCD15/{product}_2004.tif"
    epsg = get_epsg_code(sample_url)

    def do_one_year(year):
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
        return ds

    with ThreadPoolExecutor(max_workers=8) as executor:
        years = [y for y in range(2002, today.year + 1)]
        arrays = list(
            tqdm(executor.map(do_one_year, years), total=len(years))
        )

    ds = xr.concat(arrays, dim="time")
    ds.attrs["epsg"] = epsg
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


def calculate_z_score(ds, curr_month_number=None):
    curr_year = dt.datetime.now().year
    if curr_month_number is None:
        curr_month_number = dt.datetime.now().month
    clim_mean, clim_std = calculate_climatology(ds)
    curr_month = (
        ds.sel({"time": slice(f"{curr_year}-01-01", f"{curr_year}-12-31")})
        .groupby("time.month")
        .mean()
    )
    t_step = {"month": curr_month_number}
    z_score = (clim_mean.sel(t_step) - curr_month.sel(t_step)) / clim_std.sel(
        t_step
    )
    return z_score


def plot_z_score(
    ds,
    cmap=plt.cm.RdBu,
    contour=True,
    vmin=None,
    vmax=None,
    levels=np.linspace(-2.5, 2.5, 19),
):
    # The original parameters were changed to 3 degrees east and 11 degrees north
    proj = ccrs.LambertAzimuthalEqualArea(
        central_latitude=11, central_longitude=3
    )
    sns.set_context("paper")
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(8, 12))
    ax = plt.axes(projection=proj)
    # create colorbar
    cmap = plt.cm.get_cmap(cmap)
    if contour:
        norm = colors.BoundaryNorm(levels, cmap.N)
        # plot either contour of pcolormesh map
        im = ax.contourf(
            ds.coords["x"],
            ds.coords["y"],
            ds.data,
            levels=levels,
            norm=norm,
            transform=ccrs.PlateCarree(),
            spacing="uniform",
            extend="max",
            cmap=cmap,
        )
    else:
        im = ax.pcolormesh(
            ds.coords["x"],
            ds.coords["y"],
            ds.data,
            transform=ccrs.PlateCarree(),
            # norm=norm,
            # vmin=vmin, vmax=vmax,
            cmap=cmap,
        )
    ax.coastlines(resolution="10m")
    ax.add_feature(cfeature.STATES, edgecolor="gray", alpha=0.3)
    ax.add_feature(cfeature.BORDERS)
    lake = cfeature.NaturalEarthFeature(
        category="physical",
        name="lakes",
        scale="10m",
        facecolor="none",
        edgecolor="black",
    )
    ax.add_feature(lake)
    divider = make_axes_locatable(ax)
    cax = divider.append_axes(
        "right", size="5%", pad=0.05, axes_class=plt.Axes
    )
    cbar = fig.colorbar(im, cax=cax, extend="max")
    cbar.set_label("anomaly", fontsize=12)

    lon_formatter = LongitudeFormatter(zero_direction_label=False)
    ax.xaxis.set_major_formatter(lon_formatter)
    lat_formatter = LatitudeFormatter()
    ax.yaxis.set_major_formatter(lat_formatter)
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color="gray",
        alpha=1.0,
        linestyle="--",
    )
    gl.top_labels = False
    gl.right_labels = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    ax.set_extent([-3.5, 1.25, 4.5, max(ds.coords["y"])], ccrs.PlateCarree())
    
    return fig
