#!/usr/bin/env python
"""[summary]
"""
import struct
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


def write_cabo_file(site_name, lon, lat, year, c1=-0.18, c2=-0.55):
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
            "wdsp",
            "precip",
        ]
    ]
    return variables
