#!/usr/bin/env python
"""Some remote data access routines to access the NCEO
data storage on JASMIN."""
import calendar
import datetime as dt
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

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


def get_all_years(
    product,
    variable,
    first_year=2002,
    last_year=None,
    n_workers=8,
    remote_url=JASMIN_URL,
):

    if last_year is None:
        last_year = dt.datetime.today().year
    urls = {
        "MODIS": f"/vsicurl/{remote_url}/MCD15/{variable}_{first_year}wgs84.tif",
        "TAMSAT": f"/vsicurl/{remote_url}/soil_moisture/"
        + f"nc/GTiff/tamsat_{variable}_{first_year}.tif",
        "ERA5": f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_{first_year}.tif",
    }
    epsg = get_epsg_code(urls[product])

    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        years = [y for y in range(first_year, last_year + 1)]
        arrays = list(
            tqdm(
                executor.map(
                    lambda year: get_one_year(product, variable, year), years
                ),
                total=len(years),
            )
        )
    ds = xr.concat(arrays, dim="time")
    ds.attrs["epsg"] = epsg
    return ds


def get_one_year(product, variable, year, remote_url=JASMIN_URL):
    urls = {
        "MODIS": f"/vsicurl/{remote_url}/MCD15/{variable}_{year}wgs84.tif",
        "TAMSAT": f"/vsicurl/{remote_url}/soil_moisture/"
        + f"nc/GTiff/tamsat_{variable}_{year}.tif",
        "ERA5": f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_{year}.tif",
    }
    url = urls[product]
    g = gdal.Open(url, gdal.GA_ReadOnly)

    n_bands = g.RasterCount

    def to_dates(meta):
        return pd.to_datetime(f'{year}/{meta["DoY"]}', format="%Y/%j")

    def get_meta(band):
        return g.GetRasterBand(band).GetMetadata()

    dates = [to_dates(get_meta(i)) for i in range(1, n_bands + 1)]

    ds = xr.open_rasterio(url, chunks={"x": 64, "y": 64})
    ds = ds.rename({"band": "time"})
    ds = ds.assign_coords({"time": dates})
    # Scale MODIS products
    if variable == "Fpar_500m":
        ds = xr.where(ds <= 100, ds / 100.0, np.nan)
    elif variable == "Lai_500m":
        ds = xr.where(ds <= 10, ds / 100.0, np.nan)
    return ds


def check_dates(remote_url=JASMIN_URL):
    year = dt.datetime.now().year
    urls = {
        "MODIS": f"/vsicurl/{remote_url}/MCD15/Fpar_500m_{year}wgs84.tif",
        "TAMSAT": f"/vsicurl/{remote_url}/soil_moisture/"
        + f"nc/GTiff/tamsat_runoff_{year}.tif",
        "ERA5": f"/vsicurl/{remote_url}/ERA5_meteo/ssrd_{year}.tif",
    }
    dates = {}
    for product, url in urls.items():
        g = gdal.Open(url, gdal.GA_ReadOnly)
        n_bands = g.RasterCount
        meta = g.GetRasterBand(n_bands).GetMetadata()
        the_date = pd.to_datetime(f'{year}/{meta["DoY"]}', format="%Y/%j")
        the_last_month = the_date.month
        if the_date.day < calendar.monthrange(year, the_last_month)[1]:
            the_last_month -= 1

        dates[product] = the_last_month
    return dates


def add_logo(
    logo="gssti_nceo_logo2.png", origin="upper", x_o=0, y_o=0, alpha=0.5
):
    """
    Function that adds GSSTI and NCEO logo to a figure
    :param fig: Figure object to add logo to
    :param x_o: xo position of logo on figure (float)
    :param y_o: yo position of logo on figure (float)
    :return: 'logo added' (str)
    """
    logo_loc = (Path().cwd()) / logo
    if not logo_loc.exists():
        logo = [f for f in Path().cwd().rglob(f"**/{logo}")]
        logo = logo[0]
    else:
        logo = logo_loc.as_posix()

    ax = plt.axes(
        [0.6, 0.05, 0.4, 0.075], frameon=True
    )  # Change the numbers in this array to position your image [left, bottom, width, height])
    im = ax.imshow(mpimg.imread(logo))
    im.set_zorder(0)
    ax.axis("off")  # get rid of the ticks and ticklabels
    # fig.figimage(mpimg.imread(logo),
    #             xo=x_o, yo=y_o)


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
    if product.upper() == "ERA5":
        product = "ERA"
    if not product.upper() in ["ERA", "TAMSAT", "MODIS"]:
        raise ValueError(
            f'{product} isn\'t one of {["ERA", "TAMSAT", "MODIS"]}'
        )
    if product == "ERA":
        if variable not in ERA5_VARIABLES:
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
        if variable not in TAMSAT_VARIABLES:
            raise ValueError(
                f"TAMSAT product only has variables {TAMSAT_VARIABLES}"
            )
        mean_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_mean_{variable}_{period}.tif"
        std_url = f"/vsicurl/{url}/soil_moisture/nc/GTiff/clim_std_{variable}_{period}.tif"
        mean = xr.open_rasterio(mean_url, chunks={"band": 1})
        std = xr.open_rasterio(std_url, chunks={"band": 1})

    elif product == "MODIS":
        if variable not in MODIS_VARIABLES:
            raise ValueError(
                f"MODIS product only has variables {MODIS_VARIABLES}"
            )
        mean_url = f"/vsicurl/{url}/MCD15/clim_mean_{variable}_{period}.tif"
        std_url = f"/vsicurl/{url}/MCD15/clim_std_{variable}_{period}.tif"
        mean = xr.open_rasterio(
            mean_url, chunks={"band": 1, "x": 256, "y": 256}
        )
        std = xr.open_rasterio(
            std_url, chunks={"band": 1, "x": 256, "y": 256}
        )
    mean = mean.rename({"band": "month"})
    std = std.rename({"band": "month"})

    return mean, std


def calculate_climatology(
    ds, variable, first_year=2002, period="time.month", last_year=None
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
    if variable in ["precip", "rfe_filled", "runoff"]:
        monthly_ds = (
            ds.sel(time=slice(f"{first_year}-01-01", f"{last_year}-01-01"))
            .resample({"time": "1MS"})
            .sum()
        )
    else:
        monthly_ds = (
            ds.sel(time=slice(f"{first_year}-01-01", f"{last_year}-01-01"))
            .resample({"time": "1MS"})
            .mean()
        )
    clim_mean = monthly_ds.groupby(period).mean("time")
    clim_std = monthly_ds.groupby(period).std("time")
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


def do_map(
    field,
    contour=True,
    with_logo=True,
    cmap="cividis",
    vmin=-3,
    vmax=3,
    n_levels=19,
):
    sns.set_context("paper")
    sns.set_style("whitegrid")
    fig = plt.figure(figsize=(8, 12))
    cmap = plt.get_cmap(cmap)
    levels = np.linspace(vmin, vmax, n_levels)
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    if contour:
        norm = colors.BoundaryNorm(levels, cmap.N)
        # plot either contour of pcolormesh map
        im = ax.contourf(
            field.coords["x"],
            field.coords["y"],
            field.data,
            levels=levels,
            norm=norm,
            transform=ccrs.PlateCarree(),
            spacing="uniform",
            extend="max",
            cmap=cmap,
        )
    else:
        im = ax.pcolormesh(
            field.coords["x"],
            field.coords["y"],
            field.data,
            transform=ccrs.PlateCarree(),
            #            norm=norm,
            vmin=vmin,
            vmax=vmax,
            cmap=cmap,
        )
    ax.gridlines()
    ax.coastlines(resolution="50m")
    ax.add_feature(cfeature.STATES, edgecolor="gray", alpha=0.3)
    ax.add_feature(cfeature.BORDERS)
    lake = cfeature.NaturalEarthFeature(
        category="physical",
        name="lakes",
        scale="50m",
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

    ax.set_extent([-3.5, 1.25, 4.5, 11.7], ccrs.PlateCarree())
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
    if with_logo:
        add_logo()
    return fig
