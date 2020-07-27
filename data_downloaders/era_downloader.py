#!/usr/bin/env python
"""Some convenience files to grab meteorological data and convert
to CABO format to use within WOFOST. So far, using ERA5
"""
import datetime as dt
import logging
import optparse
from collections import namedtuple
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

import cdsapi
import numpy as np

JASMIN_ERA5 = "/gws/nopw/j04/odanceo/public/ERA5_meteo/"

ERAPARAMS = namedtuple(
    "ERAPARAMS", ["ssrd", "mx2t", "mn2t", "tp", "u10", "v10", "d2m"]
)
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

HELP_TEXT = """
SYNOPSIS
./era_downloader.py 
DESCRIPTION
A program to download Copernicus data.
EXAMPLES
./era_downloader.py  -v \
          

EXIT STATUS
    No exit status yet, can't be bothered.
AUTHOR
    J Gomez-Dans <j.gomez-dans@ucl.ac.uk>
    See also https://github.com/NCEO-ODA/ghana_data
"""


def grab_era5(month, year, output_folder, region, mylat, mylon):
    """A function to download ERA5 Lande data for one month. Assumes
    the Copernicus Data Service API works and is properly configured.

    Downloaded files have names `{region}.{year}_{month}.nc`.

    Function checks whether the file already exists before requesting
    the data, as requesting the data takes a while.
    
    Parameters
    ----------
    month: int
        The month number
    year: int
        The year number
    output_folder : str
        The output folder where the files will be saved to.
    region: str
        Name of the region. Useful for saving the data.
    mylat : 2-tuple, 2-list
        North and South Latitudes in decimal degrees.
    mylon : float
        West and East Longitudes in decimal degrees.
    """

    output_folder = Path(output_folder) / "netcdf"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_nc_file = (
        output_folder / f"ERA5_{region:s}.{year:d}_{month:02d}.nc"
    )
    # This is needed to keep getting the updated ERA5 datasets
    today = dt.datetime.now()
    delta_t = today - dt.datetime(year, month, 1)

    if not output_nc_file.exists() or (0 <= delta_t.days <= 120):

        LOG.info(f"Downloading {year}-{month}")
        #'80/-50/-25/0', # North, West, South, East.
        area = (
            f"{int(mylat[1]):d}/{int(mylon[0]):d}/"
            + f"{int(mylat[0]):d}/{int(mylon[1]):d}"
        )
        c = cdsapi.Client()
        c.retrieve(
            "reanalysis-era5-single-levels",
            {
                "format": "netcdf",
                "variable": [
                    "10m_u_component_of_wind",
                    "10m_v_component_of_wind",
                    "2m_dewpoint_temperature",
                    "2m_temperature",
                    "evaporation",
                    "potential_evaporation",
                    "surface_solar_radiation_downwards",
                    "total_precipitation",
                    "volumetric_soil_water_layer_1",
                    "volumetric_soil_water_layer_2",
                    "volumetric_soil_water_layer_3",
                    "volumetric_soil_water_layer_4",
                ],
                "product_type": "reanalysis",
                "year": f"{year:d}",
                "month": f"{month:02d}",
                "day": [f"{day:02d}" for day in range(1, 32)],
                "time": [f"{hour:02d}:00" for hour in range(0, 24)],
                "area": area,
                "format": "netcdf",
            },
            output_nc_file.as_posix(),
        )
        return output_nc_file.as_posix()
    else:

        LOG.info(f"Skipping {year}-{month}")
        return None


def main():
    parser = optparse.OptionParser(
        formatter=optparse.TitledHelpFormatter(), usage=HELP_TEXT
    )
    parser.add_option(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
        help="verbose output",
    )
    parser.add_option(
        "-d",
        "--data_folder",
        action="store",
        default=JASMIN_ERA5,
        dest="output_folder",
        help="Output folder to save data",
    )
    parser.add_option(
        "-r",
        "--region",
        action="store",
        dest="region",
        default="Ghana",
        help="Region name",
    )
    parser.add_option(
        "-y",
        "--lat",
        action="store",
        default="1,12",
        type=str,
        help="Minimum/maximum latitude in decimal degrees.",
    )
    parser.add_option(
        "-x",
        "--lon",
        action="store",
        default="-4,5",
        type=str,
        help="Minimum/maximum longitude in decimal degrees.",
    )

    (options, args) = parser.parse_args()
    if options.verbose:
        LOG.setLevel(logging.DEBUG)
    else:
        LOG.setLevel(logging.INFO)

    lats = [float(x) for x in options.lat.split(",")]
    lons = [float(x) for x in options.lon.split(",")]
    min_lat = min(lats)  # south
    max_lat = max(lats)  # north
    min_lon = min(lons)  # west
    max_lon = max(lons)  # east
    start_era_download(
        options.region,
        options.output_folder,
        min_lon,
        max_lon,
        min_lat,
        max_lat,
    )


def start_era_download(
    region, output_folder, min_lon, max_lon, min_lat, max_lat
):
    LOG.debug(f"Region: {region}")
    LOG.debug(f"Output: {output_folder}")
    LOG.debug(f"Lon: {min_lon}/{max_lon}")
    LOG.debug(f"Lat: {min_lat}/{max_lat}")

    wrapper = partial(
        grab_era5,
        region=region,
        output_folder=output_folder,
        mylat=[min_lat, max_lat],
        mylon=[min_lon, max_lon],
    )

    # create a thread pool of 8 threads
    years = np.arange(2000, dt.datetime.now().year + 1).astype(np.int)
    months = np.arange(1, 13).astype(np.int)
    with ThreadPoolExecutor(max_workers=8) as executor:
        for year in years:
            for month in months:
                executor.submit(wrapper, month, year)


if __name__ == "__main__":
    main()
