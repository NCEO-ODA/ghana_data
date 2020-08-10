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
import matplotlib.image as mpimg
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
TAMSAT_VARIABLES = [
    "ecan_gb",
    "esoil_gb",
    "precip",
    "runoff",
    "smc_avail_top",
    "smcl_1",
    "smcl_2",
    "smcl_3",
    "smcl_4",
]


def add_logo(fig, logo="gssti_nceo_logo2.png", x_o=400, y_o=25):
    """
    Function that adds GSSTI and NCEO logo to a figure
    :param fig: Figure object to add logo to
    :param x_o: xo position of logo on figure (float)
    :param y_o: yo position of logo on figure (float)
    :return: 'logo added' (str)
    """
    logo_arr = mpimg.imread(logo)
    fig.figimage(logo_arr, xo=x_o, yo=y_o)


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


def get_climatology(product, variable, url=JASMIN_URL, period="long"):
    if not product.upper() in ["ERA", "TAMSAT", "MODIS"]:
        raise ValueError(
            f'{product} isn\'t one of {["ERA", "TAMSAT", "MODIS"]}'
        )
    if product == "ERA":
        if not variable in ERA5_VARIABLES:
            raise ValueError(
                f"ERA5 product only has variables {ERA5_VARIABLES}"
            )
        mean_url = (
            f"/vsicurl/{url}/ERA5_meteo/clim_mean_{variable}_{period}.tif"
        )
        std_url = (
            f"/vsicurl/{url}/ERA5_meteo/clim_std_{variable}_{period}.tif"
        )
        mean = xr.open_rasterio(mean_url, chunks={"band": 1})
        std = xr.open_rasterio(std_url, chunks={"band": 1})

    elif product == "TAMSAT":
        if not variable in TAMSAT_VARIABLES:
            raise ValueError(
                f"TAMSAT product only has variables {TAMSAT_VARIABLES}"
            )
        mean_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_mean_{variable}_{period}.tif"
        std_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_std_{variable}_{period}.tif"
        mean = xr.open_rasterio(mean_url, chunks={"band": 1})
        std = xr.open_rasterio(std_url, chunks={"band": 1})

    elif product == "MODIS":
        if not variable in MODIS_VARIABLES:
            raise ValueError(
                f"MODIS product only has variables {MODIS_VARIABLES}"
            )
        mean_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_mean_{variable}_{period}.tif"
        std_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_std_{variable}_{period}.tif"
        mean = xr.open_rasterio(
            mean_url, chunks={"band": 1, "x": 256, "y": 256}
        )
        std = xr.open_rasterio(
            std_url, chunks={"band": 1, "x": 256, "y": 256}
        )
    mean = mean.rename({"band": "month"})
    std = std.rename({"band": "month"})

    return mean, std


def get_tamsat_ds(variable, remote_url=JASMIN_URL):
    assert (
        variable in TAMSAT_VARIABLES
    ), f"{variable} not one of {TAMSAT_VARIABLES}"
    today = dt.datetime.now()

    arrays = []
    sample_url = f"/vsicurl/{remote_url}/soil_moisture/nc/GTiff/tamsat_{variable}_2004.tif"
    epsg = get_epsg_code(sample_url)

    def do_one_year(year):
        url = (
            f"/vsicurl/{remote_url}/soil_moisture/"
            + f"nc/GTiff/tamsat_{variable}_{year}.tif"
        )
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


def get_modis_ds(remote_url=JASMIN_URL, product="Fpar_500m", n_workers=8):
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
        ds = xr.open_rasterio(url, chunks={"band": 1, "x": 256, "y": 256})
        ds = ds.rename({"band": "time"})
        ds = ds.assign_coords({"time": dates})
        return ds

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
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


def calculate_z_score(
    ds, clim_mean=None, clim_std=None, curr_month_number=None
):
    curr_year = dt.datetime.now().year
    if curr_month_number is None:
        curr_month_number = dt.datetime.now().month
    if clim_mean is None or clim_std is None:
        clim_mean, clim_std = calculate_climatology(ds)
    curr_month = (
        ds.sel({"time": slice(f"{curr_year}-01-01", f"{curr_year}-12-31")})
        .groupby("time.month")
        .mean()
    )
    max_month_pres = curr_month.coords["month"].values[-1]
    if curr_month_number is None:
        print(f"Last month in dataset is {max_month_pres}, using as current")
        curr_month_number = max_month_pres
    elif curr_month_number > max_month_pres:
        print(
            f"Only have data up to month {max_month_pres}, using as current"
        )
        curr_month_number = max_month_pres

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
    logo=True,
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
    ax.set_extent([-3.5, 1.25, 4.5, 11.7], ccrs.PlateCarree())
    if logo:
        add_logo()
    return fig
