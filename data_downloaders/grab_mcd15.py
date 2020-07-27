#!/usr/bin/env python
import datetime as dt
import logging
import os
from pathlib import Path

import gdal

import click

from modis_downloader import get_modis_data

from to_tif import do_tifs

try:
    MODIS_USERNAME = os.environ["MODIS_USERNAME"]
except KeyError:
    MODIS_USERNAME = ""
try:
    MODIS_PASSWORD = os.environ["MODIS_PASSWORD"]
except KeyError:
    MODIS_PASSWORD = ""

TILES = ["h17v07", "h18v07", "h17v08", "h18v08"]

MCD15_LOCATION = Path("/neodc/modis/data/MCD15A2H/collection6/")

PROCESS_LOCATION = Path("/gws/nopw/j04/odanceo/public/MCD15")

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


def download_nasa(last_doy, year, product="MCD15A2H.006"):
    username = MODIS_USERNAME
    password = MODIS_PASSWORD

    year, month, day = [
        int(x)
        for x in dt.datetime.strptime(f"{year}{last_doy}", "%Y%j")
        .strftime("%Y-%m-%d")
        .split("-")
    ]

    # You need four tiles to cover the entire Ghana
    tiles = TILES
    get_modis_data(
        username,
        password,
        "MOTA",
        product,
        tiles,
        (PROCESS_LOCATION / "hdfs").as_posix(),
        dt.datetime(year, month, day),
        n_threads=3,
    )

    logging.getLogger().setLevel(logging.INFO)


def link_neodc_files(
    curr_year,
    process_location=PROCESS_LOCATION,
    mcd_location=MCD15_LOCATION,
    product="MCD15A2H.006",
):
    """Creates symbolic links from the NEODC MODIS MCD15
    archive to the local processing archive for one year."""

    files = [
        f
        for f in (mcd_location / f"{curr_year}").rglob("MCD15A2H*.hdf")
        if f.name.split(".")[2] in TILES
    ]

    for fich in files:
        target = process_location / "hdfs" / f"{fich.name}"
        if not target.exists():
            target.symlink_to(fich)


def scan_current_files(data_loc, curr_year):
    data_loc = Path(data_loc)
    last_tif = data_loc / f"Lai_500m_{curr_year}.tif"
    g = gdal.Open(last_tif.as_posix())
    last_doy = int(g.GetRasterBand(g.RasterCount).GetMetadata()["DoY"])
    if last_doy >= 361:
        LOG.info("All done, nowt to do")
        return None
    else:
        last_time = last_doy + 1
        return last_time


@click.command()
@click.option('--location', default=PROCESS_LOCATION,
             help='Where will MODIS data be saved to')
@click.option('--product', default="MCD15A2H",
             help='MODIS product name')
@click.option('--layers', default="Lai_500m,Fpar_500m,FparLai_QC",
             help='Product layers, comma separated')
@click.option('--username', default=MODIS_USERNAME,
            help="MODIS username")
@click.option('--password', default=MODIS_PASSWORD,
            help="MODIS password")
def main(location, product, layers, username, password):
    import pdb;pdb.set_trace()
    if username == "":
        raise ValueError("No NASA username set!")
    if password == "":
        raise ValueError("No NASA username set!")
    layers = layers.split(",")
    today = dt.datetime.now()
    LOG.info(f"Started running. Current DoY {today.strftime('%Y-%j')}")
    current_year = today.year
    # Copy any files already available in JASMIN
    link_neodc_files(current_year)
    # Scan local files to see what's the latest we've processed
    last_time = scan_current_files(PROCESS_LOCATION, current_year)
    LOG.info(f"Last DoY: {last_time}")
    if dt.datetime.strptime(f"{current_year}/{last_time}", "%Y/%j") <= today:
        download_nasa(last_time, current_year)
    # Scan local files to see what's the latest we've processed
    last_time = scan_current_files(PROCESS_LOCATION, current_year)
    do_tifs(current_year, last_time, folder=PROCESS_LOCATION,
            product=product,
            layers=layers
            )

if __name__ == "__main__":
    main()