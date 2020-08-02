import subprocess

import gdal
import numpy as np
import pandas as pd
import xarray as xr

gdal.UseExceptions()

geoT = [-17.875, 0.25, 0, 37.375, 0.0, -0.25]

proj = 'GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]]'


for year in range(2000, 2021):
    cmd = [
        "cdo",
        "cat",
        f"/gws/nopw/j04/odanceo/public/soil_moisture/nc/tamsat_nceo_da_wb_v1p0p0_{year}????.nc",
        f"/gws/nopw/j04/odanceo/public/soil_moisture/nc/{year}.nc",
    ]
    with subprocess.Popen(
        cmd, stdout=subprocess.PIPE, bufsize=1, universal_newlines=True
    ) as p:
        for line in p.stdout:
            print(line, end="")  # process line here
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, p.args)

    ds = xr.open_dataset(f"{year}.nc", chunks={})

    tsteps = [
        x.strftime("%Y-%m-%d")
        for x in pd.to_datetime(ds.coords["time"].values)
    ]
    doys = [
        int(x.strftime("%j"))
        for x in pd.to_datetime(ds.coords["time"].values)
    ]

    for variable in [
        "precip",
        "esoil_gb",
        "ecan_gb",
        "smc_avail_top",
        "runoff",
    ]:
        drv = gdal.GetDriverByName("GTiff")
        dst_ds = drv.Create(
            f"/gws/nopw/j04/odanceo/public/sm_tiff/tamsat_{variable}_{year}.tif",
            278,
            292,
            len(doys),
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

        for band, (tstep, doy) in enumerate(zip(tsteps, doys)):
            this_band = dst_ds.GetRasterBand(band + 1)
            this_band.SetMetadata({"DoY": "%03d" % doy, "Date": tstep})
            x = ds[variable].isel({"time": band}).values
            x[np.isnan(x)] = -9999
            this_band.WriteArray(x)
            this_band.SetNoDataValue(-9999)
        dst_ds = None
