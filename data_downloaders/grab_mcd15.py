#!/usr/bin/env python
import datetime as dt

from pathlib import Path 

import gdal

from modis_downloader import get_modis_data

import logging

GHANA_TILES = ["h17v07", "h18v07", "h17v08", "h18v08"]

MCD15_LOCATION = Path("/neodc/modis/data/MCD15A2H/collection6/")

PROCESS_LOCATION = Path("/gws/nopw/j04/odanceo/public/MCD15")


def download_nasa(last_doy, year):
    username = "gomezdansj"
    password = "GeogG1222016"

    bothered = True
    year, month, day = [int(x)
                         for x in dt.datetime.strptime(f"{year}{last_doy}", 
                                                       "%Y%j").strftime(
                                                       "%Y-%m-%d").split("-")]

    product = "MCD15A2H.006"
    # You need four tiles to cover the entire Ghana
    tiles = ["h17v07", "h18v08", "h17v08", "h18v07"]
    get_modis_data(
        username,
        password,
        "MOTA",
        product,
        tiles,
        (PROCESS_LOCATION/"hdfs").as_posix(),
        dt.datetime(year, month, day),
        n_threads=3,
    )

    logging.getLogger().setLevel(logging.INFO) 



def link_neodc_files(curr_year,
                    process_location=PROCESS_LOCATION,
                    mcd_location=MCD15_LOCATION):
    """Creates symbolic links from the NEODC MODIS MCD15
    archive to the local processing archive for one year."""
    
    files =[f for f in  
            (mcd_location/f"{curr_year}").rglob("MCD15A2H*.hdf")
            if f.name.split(".")[2] in GHANA_TILES]

    for fich in files:
        target = process_location/"hdfs"/f"{fich.name}"
        if not target.exists():
            target.symlink_to(fich)
    

def scan_current_files(data_loc, curr_year):
    data_loc = Path(data_loc)
    last_tif = data_loc/f"Lai_500m_{curr_year}.tif"
    g = gdal.Open(last_tif.as_posix())
    last_doy = int(g.GetRasterBand(g.RasterCount)
                   .GetMetadata()['DoY'])
    if last_doy >= 361:
        print("All done, nowt to do")
    else:
        last_time = last_doy + 1
        
