# Downloading data

## Pre-requisites

Some of the services that we'll leech the data from require annoying accounts. These have to be set-up first and are assumed to be working.

### Copernicus Data Service API

Basically, you just need [to follow the instructions](https://cds.climate.copernicus.eu/api-how-to). You sign-up to the [Climate Data Store](https://cds.climate.copernicus.eu/), create an account, and then just follow the instructions.

### MODIS Downloading

For this, you need to create an [EarthData Login account](https://urs.earthdata.nasa.gov/users/new). Once you have this, you probably want to try to download a product with that account to activate it. You can do this by just clicking on this **link** and downloading it.

## Downloading the data itself

Once you have set-up the accounts above, you can easily download the data using command line tools, or developing your own codes. Bear in mind that the amount of data can be significant, so a large enough storage space is needed, as well as a reasonable Internet connection.

## Downloading the ERA5 data

The ERA5 data is the new generation ECMWF reanalysis product. As a reanalysis product, it is geared towards the past (providing an optimal estimate of the state using all available observations), and thus its usage in NRT is actually supplemented by using data from other ECMWF archives. The variables selected are:
* `10m_u_component_of_wind`
* `10m_v_component_of_wind`
* `2m_dewpoint_temperature`
* `2m_temperature`
* `evaporation`
* `potential_evaporation`
* `surface_solar_radiation_downwards`
* `total_precipitation`
* `volumetric_soil_water_layer_1`
* `volumetric_soil_water_layer_2`
* `volumetric_soil_water_layer_3`
* `volumetric_soil_water_layer_4`

(other variables can be selected by changing the code).Once the the CDS API account has been created, you can just use the simple script `data_downloaders/era_downloader.py`. The instructions are basically

```bash
Usage
=====
  
SYNOPSIS
./era_downloader.py 
DESCRIPTION
A program to download Copernicus data.
EXAMPLES
./era_downloader.py  -v           

EXIT STATUS
    No exit status yet, can't be bothered.
AUTHOR
    J Gomez-Dans <j.gomez-dans@ucl.ac.uk>
    See also https://github.com/NCEO-ODA/ghana_data


Options
=======
--help, -h              show this help message and exit
--verbose, -v           verbose output
--data_folder=OUTPUT_FOLDER, -d OUTPUT_FOLDER
                        Output folder to save data
--region=REGION, -r REGION
                        Region name
--lat=LAT, -y LAT       Minimum/maximum latitude in decimal degrees.
--lon=LON, -x LON       Minimum/maximum longitude in decimal degrees.

```

The script will download all the data available since 2000, but will check for existing files and will only download new data. So, to download data for Ghana (extent from -4 to 5 degrees longitude, 1 to 12 degrees latitude), we may use

```bash
    python ./era_downloader.py -v d era5_download/ -r Ghana -y 1,12 -x -4,5
```

Since this is likely to take a while, you may want to use nohup and log the output to a file:

```bash
    nohup python ./era_downloader.py -v d era5_download/ -r Ghana -y 1,12 -x -4,5 &>era_dload.log&
```

This downloads all the files (with names like `ERA5_Ghana.2002_12.nc`, where `Ghana` is the region name option, and we have year `2002` and month `12`). The files are stored under `era5_download/netcdf`, in the original NetCDF format. 

The files contain *hourly* data at the native resolution for the variables listed above. We usually require *daily data*, and we move to GeoTIFF format (a more network friendly data format), by using the `data_downloader/era_to_tif.py` script. This script basically reuses the location from the previous one (`era5_download`) and works on an annual basis:

```bash
python ./era_to_tif.py era5_download 2004
```