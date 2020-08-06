#!/usr/bin/env python

import numpy as np
import gdal

from calc_climatologies import ERA5_VARIABLES, MODIS_VARIABLES
from calc_climatologies import TAMSAT_VARIABLES

from calc_climatologies import get_modis_ds
from calc_climatologies import get_era5_ds
from calc_climatologies import get_tamsat_ds


if __name__ == "__main__":
    
    clim_periods = {'long':[2002, 2019],
                    'recent':[2015, 2019],
                    'past': [2002, 2010]
                    }
    for variable in ERA5_VARIABLES:
        ds = get_era5_ds(variable)
        for k, v in clim_periods.items():
            m, s = calculate_climatology(ds, first_year=v[0], period="time.month", last_year=v[1])
            break
        break
        
    
