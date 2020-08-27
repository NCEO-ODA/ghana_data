#!/usr/bin/env python
"""[summary]
"""
import struct
from pathlib import Path
from textwrap import dedent

import gdal
import numpy
from basic_calcs import ERA5_VARIABLES, get_all_years, get_one_year

DEM_FILE = "/vsicurl/http://www2.geog.ucl.ac.uk/~ucfafyi/eles/global_dem.vrt"


def retrieve_pixel_value(lon, lat, data_source):
    """Retrieve pixel value from a GDAL-friendly dataset.

    We assume the data type of the raster here!!!!

    Parameters
    ----------
    lon : float
        Longitude in decimal degrees
    lat : float
        Latitude in decimal degrees
    data_source : str
        An existing GDAL-readable dataset. Can be remote.

    Returns
    -------
    int
       The value of the pixel.
    """
    dataset = gdal.Open(data_source)

    gt = dataset.GetGeoTransform()
    the_band = dataset.GetRasterBand(1)
    px = int((lon - gt[0]) / gt[1])  # x pixel
    py = int((lat - gt[3]) / gt[5])  # y pixel

    buf = the_band.ReadRaster(px, py, 1, 1, buf_type=gdal.GDT_Int16)
    elev = struct.unpack("h", buf)

    return elev[0]


def get_era5_data(lon, lat, year):
    data_stacks = {
        variable: get_one_year("ERA5", variable, year)
        for variable in ERA5_VARIABLES
    }
    return data_stacks


def write_cabo_file(
    site_name, lon, lat, year, dest_folder=".", c1=-0.18, c2=-0.55
):
    cabo_file = Path(dest_folder) / f"{site_name}_{year}.txt"
    elev = retrieve_pixel_value(lon, lat, DEM_FILE)
    hdr_chunk = f"""\
    *---------------------------------------------------
    * Station: {site_name:s}
    * Year: {year:d}
    * Origin: ERA5-Reanalysis
    * Columns & units
    * ===================
    * 1. station number
    * 2. year
    * 3. Day of Year
    * 4. Irradiance   (kJ·m-2·d-1)
    * 5. Daily minimum temperature (degC)
    * 6. Daily maximum temperature (degC)
    * 7. Vapour pressure (kPa)
    * 8. Mean wind speed (m·s-1)
    * 9. Precipitation (mm·d-1)
    ** WCCDESCRIPTION={site_name:s}
    ** WCCFORMAT=2
    ** WCCYEARNR={year:d}
    *------------------------------------------------------------*
    {lon:.2f}  {lat:.2f}  {elev:.2f} {c1:.2f}  {c2:.2f}
    """
    hdr_chunk = dedent(hdr_chunk)
    data_stack = get_era5_data(lon, lat, year)
    variables = [
        data_stack[variable].sel(y=lat, x=lon, method="nearest").values
        for variable in [
            "ssrd",
            "t2m_min",
            "t2m_max",
            "hum",
            "wspd",
            "precip",
        ]
    ]
    station_number = 1
    with cabo_file.open("w") as fp:
        fp.write(hdr_chunk)
        for d in range(data_stack["ssrd"].shape[0]):
            fp.write(
                f"{station_number:d}\t{year:d}\t{d+1:d}\t"
                + f"{round(variables[0][d]):5.1f}\t"
                + f"{round(variables[1][d]*10/10):5.1f}\t"
                + f"{round(variables[2][d]*10/10):5.1f}\t"
                + f"{round(variables[3][d]*1000/1000):5.3f}\t"
                + f"{round(variables[4][d]*10/10):4.1f}\t"
                + f"{round(variables[5][d]*10/10):4.1f}\n"
            )
