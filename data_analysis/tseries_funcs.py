import datetime as dt

import gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr

from .basic_calcs import (
    ERA5_VARIABLES,
    MODIS_VARIABLES,
    TAMSAT_VARIABLES,
    calculate_z_score,
    check_dates,
    do_map,
    get_climatology,
    get_climatology_url,
    get_one_year,
)
from .simple_gui import PRODUCT_UNITS

JASMIN_URL = "http://gws-access.ceda.ac.uk/public/odanceo/"


## UNITS for different products.
## Use a string for the units, can use latex if
## needed. Keep dictionary keys the same as ERA5_VARIABLES eg

ERA5_UNITS = {
    "hum": "Vapour Pressure \n" + r"$[kPa]$",
    "precip": "Daily precipitation \n" + r"$mm\cdot d^{-1}$",
    "ssrd": "Daily irradiance \n" + r"$W\cdot m^{-2}$",
    "t2m_max": "Maximum daily temperature \n" + r"$^{\circ}C$",
    "t2m_mean": "Mean daily temperature \n" + r" $^{\circ}C$",
    "t2m_min": "Minimum daily temperature \n" + r" $^{\circ}C$",
    "wspd": "Windspeed \n" + r" $m\cdot s^{-1}$",
}

MODIS_UNITS = {
    "Fpar_500m": "$fAPAR$ \n" + r" [-]",
    "Lai_500m": "Leaf Area Index \n" + r" $m^{2}\cdot m^{-2}$",
    "FparLai_QC": "N/A",
}

TAMSAT_UNITS = {
    "ecan_gb": "Evaporation from surface store \n"
    + r"(kg m$^{-2}$ s$^{-1}$)",
    "esoil_gb": "Evaporationtranspiration \n" + r"(kg m$^{-2}$ s$^{-1}$)",
    "precip": "Precipitation rate \n"
    + r"(kg m$^{-2}$ s$^{-1}$)",  # multiply this by 86400 to get (mm day-1)
    "runoff": "Gridbox runoff rate \n" + r"(kg m$^{-2}$ s$^{-1}$)",
    "smc_avail_top": "Available soil moisture in surface layer \n"
    + r"(kg m$^{-2}$)",  # divide by 1000 to get (m3 m-3)
    "smcl_1": "Top layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 10cm thick, divide by 100 to get (m3 m-3)
    "smcl_2": "Second layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 25cm thick, divide by 250 to get (m3 m-3)
    "smcl_3": "Third layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 65cm thick, divide by 650 to get (m3 m-3)
    "smcl_4": "Fourth layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 200cm thick, divide by 2000 to get (m3 m-3)
}


variable_lists = {
    "TAMSAT": TAMSAT_VARIABLES,
    "ERA5": ERA5_VARIABLES,
    "MODIS": MODIS_VARIABLES,
}

periods = {"2002-2019": "long", "2015-2019": "recent", "2002-2010": "past"}

landcover = {
    2015: [
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2015/E000N20/E000N20_PROBAV_LC100_global_v3.0.1_2015-base_Discrete-Classification-map_EPSG-4326.tif",
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2015/W020N20/W020N20_PROBAV_LC100_global_v3.0.1_2015-base_Discrete-Classification-map_EPSG-4326.tif",
    ],
    2016: [
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2016/E000N20/E000N20_PROBAV_LC100_global_v3.0.1_2016-conso_Discrete-Classification-map_EPSG-4326.tif",
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2016/W020N20/W020N20_PROBAV_LC100_global_v3.0.1_2016-conso_Discrete-Classification-map_EPSG-4326.tif",
    ],
    2017: [
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2017/E000N20/E000N20_PROBAV_LC100_global_v3.0.1_2017-conso_Discrete-Classification-map_EPSG-4326.tif",
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2017/W020N20/W020N20_PROBAV_LC100_global_v3.0.1_2017-conso_Discrete-Classification-map_EPSG-4326.tif",
    ],
    2018: [
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2018/E000N20/E000N20_PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif",
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2018/W020N20/W020N20_PROBAV_LC100_global_v3.0.1_2018-conso_Discrete-Classification-map_EPSG-4326.tif",
    ],
    2019: [
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2019/E000N20/E000N20_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
        "https://s3-eu-west-1.amazonaws.com/vito.landcover.global/v3.0.1/2019/W020N20/W020N20_PROBAV_LC100_global_v3.0.1_2019-nrt_Discrete-Classification-map_EPSG-4326.tif",
    ],
}


def agricultural_areas(
    landcover=landcover,
    carto_folder="/gws/nopw/j04/odanceo/public/ghana_carto/",
):

    for year, files in landcover.items():
        file_ptrs = [f"/vsicurl/{fich}" for fich in files]
        _ = gdal.Warp(
            f"{carto_folder}/croplands_{year}_LC100.tif",
            file_ptrs,
            format="GTiff",
            cutlineDSName=f"{carto_folder}/GHA_admbndp1_1m_GAUL.shp",
            cropToCutline=True,
            creationOptions=[
                "TILED=YES",
                "COMPRESS=DEFLATE",
                "PREDICTOR=1",
                "BLOCKXSIZE=512",
                "BLOCKYSIZE=512",
            ],
        )


def crop_ds(url, year, outputBounds, xRes_lc, yRes_lc):

    product = gdal.Warp(
        "",
        url,
        format="MEM",
        outputBounds=outputBounds,
        xRes=xRes_lc,
        yRes=yRes_lc,
    )

    n_bands = product.RasterCount

    def to_dates(meta):
        return pd.to_datetime(f'{year}/{meta["DoY"]}', format="%Y/%j")

    def get_meta(band):
        return product.GetRasterBand(band).GetMetadata()

    data = product.ReadAsArray()

    try:
        dates = [to_dates(get_meta(i)) for i in range(1, n_bands + 1)]
    except KeyError:
        return data
    return data, dates


def to_xarray(array, dates):
    return xr.Dataset.from_dict(
        {
            "coords": {
                "t": {"dims": "t", "data": dates, "attrs": {"units": "s"}}
            },
            "dims": ["t", "x", "y"],
            "data_vars": {
                "variable": {"dims": ["t", "y", "x"], "data": array},
            },
        }
    )


def get_one_region_landcover(
    product,
    variable,
    year,
    region_where,
    region_ds="GHA_admbndp2_1m_GAUL.shp",
    ax=None,
    period="long",
    landcover_file="croplands_2018_LC100.tif",
    remote_url=JASMIN_URL,
):
    if ax is None:
        ax = plt.gca()

    landcover_file = f"/vsicurl/{remote_url}/ghana_carto/{landcover_file}"
    region_ds = f"/vsicurl/{remote_url}/ghana_carto/{region_ds}"
    clim_mean, clim_std, chunks = get_climatology_url(
        product, variable, period=period
    )

    xx = gdal.Info(landcover_file, format="json")
    geoT = xx["geoTransform"]

    xRes_lc, yRes_lc = geoT[1], np.abs(geoT[5])
    lc = gdal.Warp(
        "",
        landcover_file,
        format="MEM",
        cutlineDSName=region_ds,
        #                   cutlineSQL=region_where,
        cutlineWhere=region_where,
        cropToCutline=True,
        xRes=xRes_lc,
        yRes=yRes_lc,
    )
    meta_lc = gdal.Info(lc, format="json")
    outputBounds = (
        meta_lc["cornerCoordinates"]["lowerLeft"]
        + meta_lc["cornerCoordinates"]["upperRight"]
    )

    urls = {
        "MODIS": f"/vsicurl/{remote_url}/MCD15/{variable}_{year}wgs84.tif",
        "TAMSAT": f"/vsicurl/{remote_url}/soil_moisture/"
        + f"nc/GTiff/tamsat_{variable}_{year}.tif",
        "ERA5": f"/vsicurl/{remote_url}/ERA5_meteo/{variable}_{year}.tif",
    }
    url = urls[product]
    product_magnitude, dates = crop_ds(
        url, year, outputBounds, xRes_lc, yRes_lc
    )
    if variable == "Fpar_500m":
        product_magnitude = np.where(
            product_magnitude <= 100, product_magnitude / 100.0, np.nan
        )
    elif variable == "Lai_500m":
        product_magnitude = np.where(
            product_magnitude <= 10, product_magnitude / 100.0, np.nan
        )

    clim_mean = crop_ds(clim_mean, year, outputBounds, xRes_lc, yRes_lc)
    clim_std = crop_ds(clim_std, year, outputBounds, xRes_lc, yRes_lc)

    lc = lc.ReadAsArray()
    lc_m = np.isin(lc, [200, 0, 50, 60, 70, 80, 90, 100], invert=True)

    product_magnitude = 1.0 * product_magnitude
    product_magnitude[:, ~lc] = np.nan
    clim_mean[:, ~lc_m] = np.nan
    clim_std[:, ~lc_m] = np.nan

    ds = to_xarray(product_magnitude, dates)
    if variable in ["precip", "rfe_filled", "runoff"]:
        monthly_ds = ds.resample(
            {"t": "1MS"}
        ).sum()  # ,  loffset=pd.Timedelta(14, 'd'))
    else:
        monthly_ds = ds.resample(
            {"t": "1MS"}
        ).mean()  # ,  loffset=pd.Timedelta(14, 'd'))

    product_magnitude = monthly_ds.variable.values
    product_magnitude[:, ~lc_m] = np.nan

    dates = monthly_ds.coords["t"].values
    clim_dates = [dt.datetime(year, i, 1) for i in range(1, 13)]

    x = np.nanmean(product_magnitude, axis=(1, 2))
    y = np.nanstd(product_magnitude, axis=(1, 2))
    line = ax.plot(dates, x, "-", label=year)
    line_color = line[0].get_c()
    ax.fill_between(dates, x - y, x + y, color=line_color, alpha=0.5)
    x_clim = np.nanmean(clim_mean, axis=(1, 2))
    y_clim = np.nanmean(clim_std, axis=(1, 2))
    line = ax.plot(clim_dates, x_clim, "-", label="LTA")
    line_color = line[0].get_c()
    ax.fill_between(
        clim_dates,
        x_clim - y_clim,
        x_clim + y_clim,
        color=line_color,
        alpha=0.5,
    )
    ax.set_ylabel(PRODUCT_UNITS[product][variable], fontsize=9)
    ax.legend(loc="upper right", fontsize=9)
