#!/usr/bin/env python
from osgeo import gdal
import numpy as np
import xarray as xr
from basic_calcs import (
    ERA5_VARIABLES,
    TAMSAT_VARIABLES,
    calculate_climatology,
    get_all_years,
)

# from dask import LocalCluster


def to_tif(x, fname, nx, ny, n_bands, geoT, proj):
    drv = gdal.GetDriverByName("GTiff")
    dst_ds = drv.Create(
        fname,
        nx,
        ny,
        n_bands,
        gdal.GDT_Float32,
        options=[
            "COMPRESS=DEFLATE",
            "TILED=YES",
            "BIGTIFF=YES",
            "PREDICTOR=1",
        ],
    )
    dst_ds.SetGeoTransform(geoT)
    dst_ds.SetProjection(proj)
    for band in range(n_bands):
        print(band)
        dst_ds.GetRasterBand(band + 1).WriteArray(
            x.isel({"month": band}).values
        )
    dst_ds = None


if __name__ == "__main__":

    clim_periods = {
        "long": [2002, 2019],
        "recent": [2015, 2019],
        "past": [2002, 2010],
    }
    ERA5 = True
    TAMSAT = True
    MODIS = True
    if ERA5:
        ####### ERA5
        g = gdal.Open(
            "/gws/nopw/j04/odanceo/public/ERA5_meteo/Ghana_precip_2001.tif"
        )
        geoT = g.GetGeoTransform()
        proj = g.GetProjectionRef()
        nx, ny = g.RasterXSize, g.RasterYSize
        for variable in ERA5_VARIABLES:
            ds = get_all_years("ERA5", variable)
            for k, v in clim_periods.items():
                m, s = calculate_climatology(
                    ds,
                    variable,
                    first_year=v[0],
                    period="time.month",
                    last_year=v[1],
                )

                fname = f"/gws/nopw/j04/odanceo/public/ERA5_meteo/clim_mean_{variable}_{k}.tif"
                to_tif(m.compute(), fname, nx, ny, 12, geoT, proj)
                fname = f"/gws/nopw/j04/odanceo/public/ERA5_meteo/clim_std_{variable}_{k}.tif"
                to_tif(s.compute(), fname, nx, ny, 12, geoT, proj)

    if TAMSAT:
        ####### TAMSAT
        g = gdal.Open(
            "/gws/nopw/j04/odanceo/public/soil_moisture/nc/GTiff/tamsat_precip_2001.tif"
        )
        geoT = g.GetGeoTransform()
        proj = g.GetProjectionRef()
        nx, ny = g.RasterXSize, g.RasterYSize
        for variable in TAMSAT_VARIABLES:
            ds = get_all_years("TAMSAT", variable)
            for k, v in clim_periods.items():
                m, s = calculate_climatology(
                    ds,
                    variable,
                    first_year=v[0],
                    period="time.month",
                    last_year=v[1],
                )

                fname = f"/gws/nopw/j04/odanceo/public/soil_moisture/nc/GTiff/clim_mean_{variable}_{k}.tif"
                to_tif(m.compute(), fname, nx, ny, 12, geoT, proj)
                fname = f"/gws/nopw/j04/odanceo/public/soil_moisture/nc/GTiff/clim_std_{variable}_{k}.tif"
                to_tif(s.compute(), fname, nx, ny, 12, geoT, proj)
    if MODIS:
        ### MODIS
        g = gdal.Open(
            "/gws/nopw/j04/odanceo/public/MCD15/Lai_500m_2010wgs84.tif"
        )
        geoT = g.GetGeoTransform()
        proj = g.GetProjectionRef()
        nx, ny = g.RasterXSize, g.RasterYSize
        for variable in ["Lai_500m", "Fpar_500m"]:
            ds = get_all_years("MODIS", variable)
            for k, v in clim_periods.items():
                print(variable, k)
                m, s = calculate_climatology(
                    ds,
                    variable,
                    first_year=v[0],
                    period="time.month",
                    last_year=v[1],
                )
                print("dumping to disk...")
                fname = f"/gws/nopw/j04/odanceo/public/MCD15/clim_mean_{variable}_{k}.tif"
                to_tif(m, fname, nx, ny, 12, geoT, proj)
                fname = f"/gws/nopw/j04/odanceo/public/MCD15/clim_std_{variable}_{k}.tif"
                to_tif(s, fname, nx, ny, 12, geoT, proj)
