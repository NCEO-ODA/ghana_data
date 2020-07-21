#!/usr/bin/env python
import datetime as dt
from pathlib import Path 

import gdal
import osr
import xarray as xr
import numpy as np


def write_tif(arr, var, year, loc, geoT, srs):
    n_bands, ny, nx = arr.shape
    fname_out = (loc/f"{var}_{year}.tif").as_posix()
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(fname_out,
                        nx, ny, n_bands, gdal.GDT_Float32,
                        options=["COMPRESS=DEFLATE", "TILED=YES", 
                                 "BIGTIFF=YES", "PREDICTOR=1"])
    ds.SetGeoTransform(geoT)
    ds.SetProjection(srs)
    this_time = dt.date(year, 1, 1)
    delta = dt.timedelta(days=1)
    for band in range(n_bands):
        data_band = ds.GetRasterBand(band+1)
        data_band.SetMetadata({"DoY":f"{band+1}",
                               "Date":this_time.strftime("%Y-%m-%d")})
        data_band.WriteArray(arr[band, :, :])
        this_time += delta
    ds = None
    print(f"{fname_out} saved")

    
def to_sensible_format(loc, year, month, ):

    p = Path(loc) 
    g = gdal.Open(f'NETCDF:"{loc}/netcdf/ERA5_Ghana.{year}_{month:02d}.nc":ssrd')
    geoT = g.GetGeoTransform()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    srs = srs.ExportToWkt()
    
    nx, ny = g.RasterXSize, g.RasterYSize
    for year in range(2000, dt.datetime.now().year + 1):
        fnames = sorted([f for f in (p).glob(f"ERA5_Ghana.{year}_??.nc")])   
        ds=xr.concat([xr.open_dataset(f, chunks={}, mask_and_scale=True ) 
                    for f in fnames], 
                    "time") 
        ssrd = ds.ssrd.resample(time="1D").mean()/1000.
        write_tif(ssrd.values, "ssrd", year, p, geoT, srs)
        tp = ds.tp.resample(time="1D").sum()*1000.
        tp= xr.where(tp >= 0.001, tp, 0)
        write_tif(tp.values, "precip", year, p, geoT, srs)
        t2m = ds.t2m.resample(time="1D").mean() - 273.15
        write_tif(t2m.values, "t2m_mean", year, p, geoT, srs)
        t2m = ds.t2m.resample(time="1D").min() - 273.15
        write_tif(t2m.values, "t2m_min", year, p, geoT, srs)
        t2m = ds.t2m.resample(time="1D").max() - 273.15
        write_tif(t2m.values, "t2m_max", year, p, geoT, srs)
        tdew = ds.d2m.resample(time="1D").mean() - 273.15
        tmp = (17.27 * tdew) / (tdew + 237.3)
        ea = 0.6108 * np.exp(tmp)
        write_tif(ea.values, "hum", year, p, geoT, srs)
        u10 = ds.u10.resample(time="1D").mean()
        v10 = ds.v10.resample(time="1D").mean()
        wspd = np.sqrt(u10**2 + v10**2)
        write_tif(wspd.values, "wspd", year, p, geoT, srs)

