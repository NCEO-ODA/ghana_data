#!/usr/bin/env python
import datetime as dt
import logging
from pathlib import Path

import click
import gdal
import numpy as np
import osr
import xarray as xr

gdal.UseExceptions()

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
if not LOG.handlers:
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - " + "%(levelname)s - %(message)s"
    )
    ch.setFormatter(formatter)
    LOG.addHandler(ch)
LOG.propagate = False


def write_tif(arr, var, year, loc, geoT, srs):
    n_bands, ny, nx = arr.shape
    fname_out = (loc / f"{var}_{year}.tif").as_posix()
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(
        fname_out,
        nx,
        ny,
        n_bands,
        gdal.GDT_Float32,
        options=[
            "COMPRESS=DEFLATE",
            "TILED=YES",
            "BIGTIFF=YES",
            "PREDICTOR=1",
            "BLOCKXSIZE=32",
            "BLOCKYSIZE=32",
        ],
    )
    ds.SetGeoTransform(geoT)
    ds.SetProjection(srs)
    this_time = dt.date(year, 1, 1)
    delta = dt.timedelta(days=1)
    for band in range(n_bands):
        data_band = ds.GetRasterBand(band + 1)
        data_band.SetMetadata(
            {"DoY": f"{band+1}", "Date": this_time.strftime("%Y-%m-%d")}
        )
        data_band.WriteArray(arr[band, :, :])
        this_time += delta
    ds = None
    LOG.info(f"{fname_out} saved")


@click.command()
@click.argument("loc")
@click.argument("year")
# @click.option("--config_file", type=click.Path())
def to_sensible_format(loc, year):
    year = int(year)
    loc = Path(loc)
    LOG.info(f"Doing {year} on {loc.as_posix()}")
    if not loc.exists():
        raise IOError(f"{loc} does not exist!")
    g = gdal.Open(
        f'NETCDF:"{loc.as_posix()}/netcdf/ERA5_Ghana.{year}_01.nc":ssrd'
    )
    geoT = g.GetGeoTransform()
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    srs = srs.ExportToWkt()
    # nx, ny = g.RasterXSize, g.RasterYSize

    fnames = sorted(
        [f for f in (loc / "netcdf").glob(f"ERA5_Ghana.{year}_??.nc")]
    )
    ds = xr.concat(
        [xr.open_dataset(f, chunks={}, mask_and_scale=True) for f in fnames],
        "time",
    )
    if "expver" in ds.coords:
        LOG.info("ERA5RT data. Expunging one dimension, hope it works")
        ds = ds.drop_sel({"expver": 1}).squeeze()
    ssrd = ds.ssrd.resample(time="1D").sum() / 1000.0
    write_tif(ssrd.values, "ssrd", year, loc, geoT, srs)
    tp = ds.tp.resample(time="1D").sum() * 1000.0
    tp = xr.where(tp >= 0.001, tp, 0)
    write_tif(tp.values, "precip", year, loc, geoT, srs)
    t2m = ds.t2m.resample(time="1D").mean() - 273.15
    write_tif(t2m.values, "t2m_mean", year, loc, geoT, srs)
    t2m = ds.t2m.resample(time="1D").min() - 273.15
    write_tif(t2m.values, "t2m_min", year, loc, geoT, srs)
    t2m = ds.t2m.resample(time="1D").max() - 273.15
    write_tif(t2m.values, "t2m_max", year, loc, geoT, srs)
    tdew = ds.d2m.resample(time="1D").mean() - 273.15
    tmp = (17.27 * tdew) / (tdew + 237.3)
    ea = 0.6108 * np.exp(tmp)
    write_tif(ea.values, "hum", year, loc, geoT, srs)
    u10 = ds.u10.resample(time="1D").mean()
    v10 = ds.v10.resample(time="1D").mean()
    wspd = np.sqrt(u10 ** 2 + v10 ** 2)
    write_tif(wspd.values, "wspd", year, loc, geoT, srs)
    LOG.info("Successfully done!")


if __name__ == "__main__":
    to_sensible_format()
