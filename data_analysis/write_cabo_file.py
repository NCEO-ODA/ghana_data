#!/usr/bin/env python
"""[summary]
"""
import struct
from pathlib import Path
from textwrap import dedent

import gdal
import numpy as np
import pandas as pd

from .basic_calcs import ERA5_VARIABLES, get_all_years, get_one_year

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


def get_era5_data(years):
    if type(years) != list:
        years = [
            years,
        ]
    if len(years) == 1:
        data_stacks = {
            variable: get_one_year("ERA5", variable, years[0])
            for variable in ERA5_VARIABLES
        }
    else:
        data_stacks = {
            variable: get_all_years(
                "ERA5", variable, first_year=min(years), last_year=max(years)
            )
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
    data_stack = get_era5_data(year)
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


def write_pcse_csv(
    country, site_name, lon, lat, years, dest_folder=".", c1=-0.18, c2=-0.55
):
    if type(years) != list:
        years = list(years)
    fyears = "_".join([f"{y}" for y in years])
    csv_file = Path(dest_folder) / f"{country}_{site_name}_{fyears}.csv"
    elev = retrieve_pixel_value(lon, lat, DEM_FILE)

    hdr_chunk = f"""Country     = '{country}'
Station     = '{site_name}'
Description = 'Reanalysis data'
Source      = 'ERA5'
Contact     = 'J Gomez-Dans'
Longitude = {lon}; Latitude = {lat}; Elevation = {elev}; AngstromA = {c1}; AngstromB = {c2}; HasSunshine = False
## Daily weather observations (missing values are NaN)
    """

    hdr_chunk = dedent(hdr_chunk)
    data_stack = get_era5_data(years)
    variables = {
        "time": data_stack["ssrd"]
        .coords["time"]
        .values.astype("M8[D]")
        .astype("O")
    }
    for variable in ["ssrd", "t2m_min", "t2m_max", "hum", "wspd", "precip"]:
        variables[variable] = (
            data_stack[variable].sel(y=lat, x=lon, method="nearest").values
        )
    variables["SNOWDEPATH"] = variables["ssrd"] * np.nan
    variables = pd.DataFrame(variables)
    variables.columns = [
        "DAY",
        "IRRAD",
        "TMIN",
        "TMAX",
        "VAP",
        "WIND",
        "RAIN",
        "SNOWDEPTH",
    ]
    with csv_file.open(mode="w", newline="") as fp:
        fp.write(hdr_chunk)
        variables.to_csv(fp, index=False, na_rep="nan")
    print(f"Saved to {csv_file}")


#    20040101,NaN,-0.7,1.1,0.55,3.6,0.5,NaN
#    20040102,3888,-7.5,0.9,0.44,3.1,0,NaN
#    20040103,2074,-6.8,-0.5,0.45,1.8,0,NaN
#    20040104,1814,-3.6,5.9,0.66,3.2,2.5,NaN
#    20040105,1469,3,5.7,0.78,2.3,1.3,NaN
#    [...]
