#!/usr/bin/env python
import shutil
import sys
from pathlib import Path

import gdal
import numpy as np


def get_sfc_qc(qa_data, mask57=0b11100000):
    sfc_qa = np.right_shift(np.bitwise_and(qa_data, mask57), 5)
    return sfc_qa


def mosaic_dates_wgs84_country(
    year,
    doy,
    folder="data",
    product="MOD13A2",
    country="Ghana",
    layer="Lai_500m",
    srs="EPSG:2136",
    fmt="MEM",
):
    """Searches for all files in folder `folder` of a given MODIS product for
    the selected year and day of year, and mosaics them together, clipping by some
    shapefile. Returns the actual raster data, as well as the corner coordinates in 
    WGS-84 projection.
    """
    folder = Path(folder)
    files = [
        f'HDF4_EOS:EOS_GRID:"{f.as_posix():s}":MOD_Grid_MOD15A2H:{layer:s}'
        for f in folder.glob(f"{product:s}.A{year:04d}{doy:03d}.*hdf")
    ]
    if files:
        g = gdal.Warp(
            (folder / f"{layer:s}.A{year:04d}{doy:03d}.tif").as_posix(),
            files,
            format=fmt,
            dstSRS=srs,
            xRes=500,
            yRes=500,
            cutlineDSName="/gws/nopw/j04/odanceo/public/MCD15/ne_50m_admin_0_countries.shp",
            cutlineWhere=f"NAME='{country:s}'",
            cropToCutline=True,
            creationOptions=[
                "TILED=YES",
                "INTERLEAVE=BAND",
                "COMPRESS=DEFLATE",
                "COPY_SRC_OVERVIEWS=YES",
            ],
        )
        g = None
        return folder / f"{layer:s}.A{year:04d}{doy:03d}.tif"
    else:
        return None


def do_tifs(year, max_doy, folder="/gws/nopw/j04/odanceo/public/MCD15"):

    folder = Path(folder)
    hdf_folder = folder / "hdfs"
    for layer in ["Lai_500m", "Fpar_500m", "FparLai_QC"]:
        fnames = []
        for doy in range(1, max_doy + 1, 8):
            fname = mosaic_dates_wgs84_country(
                year,
                doy,
                hdf_folder.as_posix(),
                product="MCD15A2H",
                layer=layer,
                fmt="GTiff",
            )
            if fname is not None:
                fnames.append(fname)

        fnames = sorted([f.as_posix() for f in fnames])
        dst_ds = gdal.BuildVRT(
            (folder / f"{layer:s}_{year:d}.vrt").as_posix(),
            fnames,
            options=gdal.BuildVRTOptions(separate=True),
        )
        dst_ds = None
        dst_ds = gdal.Translate(
            (folder / f"{layer:s}_{year:d}.tif").as_posix(),
            (folder / f"{layer:s}_{year:d}.vrt").as_posix(),
            options=gdal.TranslateOptions(
                format="GTiff",
                creationOptions=[
                    "TILED=YES",
                    "INTERLEAVE=BAND",
                    "COMPRESS=LZW",
                    "COPY_SRC_OVERVIEWS=YES",
                ],
            ),
        )
        doy = 1
        for band in range(1, dst_ds.RasterCount + 1):
            dst_ds.GetRasterBand(band).SetMetadata({"DoY": f"{doy:d}"})
            doy += 8
        dst_ds = None
        g = gdal.Open(
            (folder / f"{layer:s}_{year:d}.tif").as_posix(), gdal.GA_Update
        )
        g.BuildOverviews("average", np.power(2, np.arange(8)))
        g = None
        g = gdal.Translate(
            (folder / f"{layer:s}_{year:d}_temporary.tif").as_posix(),
            (folder / f"{layer:s}_{year:d}.tif").as_posix(),
            format="GTiff",
            creationOptions=[
                "TILED=YES",
                "INTERLEAVE=BAND",
                "COMPRESS=LZW",
                "COPY_SRC_OVERVIEWS=YES",
            ],
        )
        shutil.move(
            folder / f"{layer:s}_{year:d}_temporary.tif",
            (folder / f"{layer:s}_{year:d}.tif").as_posix(),
        )


# if __name__ == "__main__":
#    do_tifs(2019, 23)
