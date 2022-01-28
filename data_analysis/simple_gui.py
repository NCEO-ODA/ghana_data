#!/usr/bin/env python
"""A simple Jupyter notebook GUI to do spatial plots over Ghana.
"""
import datetime as dt

from .basic_calcs import (
    ERA5_VARIABLES,
    MODIS_VARIABLES,
    TAMSAT_VARIABLES,
    calculate_z_score,
    check_dates,
    do_map,
    get_climatology,
    get_one_year,
)

try:
    import ipywidgets.widgets as widgets
except ImportError:
    print("You will not be able to use the interactive GUI")


## UNITS for different products.
## Use a string for the units, can use latex if
## needed. Keep dictionary keys the same as ERA5_VARIABLES eg

ERA5_UNITS = {
    "hum": "Vapour Pressure \n" + r"$[kPa]$",
    "precip": "Daily precipitation \n" + r"$mm\cdot d^{-1}$",
    "ssrd": "Daily irradiance \n" + r"$W\cdot m^{-2}$",
    "t2m_max": "Maximum daily temperature \n" + r"$^{\circ}C$",
    "t2m_mean": "Mean daily temperature \n" + r" $^{\circ}C$",
    "t2m_min": "Minimum daily temperature \n" + r" $^{\circ}C$",
    "wspd": "Windspeed \n" + r" $m\cdot s^{-1}$",
}

MODIS_UNITS = {
    "Fpar_500m": "$fAPAR$ \n" + r" [-]",
    "Lai_500m": "Leaf Area Index \n" + r" $m^{2}\cdot m^{-2}$",
    "FparLai_QC": "N/A",
}

TAMSAT_UNITS = {
    "ecan_gb": "Evaporation from surface store \n"
    + r"(kg m$^{-2}$ s$^{-1}$)",
    "esoil_gb": "Evaporationtranspiration \n" + r"(kg m$^{-2}$ s$^{-1}$)",
    "precip": "Precipitation rate \n"
    + r"(kg m$^{-2}$ s$^{-1}$)",  # multiply this by 86400 to get (mm day-1)
    "runoff": "Gridbox runoff rate \n" + r"(kg m$^{-2}$ s$^{-1}$)",
    "smc_avail_top": "Available soil moisture in surface layer \n"
    + r"(kg m$^{-2}$)",  # divide by 1000 to get (m3 m-3)
    "smcl_1": "Top layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 10cm thick, divide by 100 to get (m3 m-3)
    "smcl_2": "Second layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 25cm thick, divide by 250 to get (m3 m-3)
    "smcl_3": "Third layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 65cm thick, divide by 650 to get (m3 m-3)
    "smcl_4": "Fourth layer soil moisture \n"
    + r"(kg m$^{-2}$)",  # 200cm thick, divide by 2000 to get (m3 m-3)
}

PRODUCT_UNITS = {
    "ERA5": ERA5_UNITS,
    "TAMSAT": TAMSAT_UNITS,
    "MODIS": MODIS_UNITS,
}

variable_lists = {
    "TAMSAT": TAMSAT_VARIABLES,
    "ERA5": ERA5_VARIABLES,
    "MODIS": MODIS_VARIABLES,
}

periods = {"2002-2019": "long", "2015-2019": "recent", "2002-2010": "past"}

diverging_cmaps = [
    "PiYG",
    "PRGn",
    "BrBG",
    "PuOr",
    "RdGy",
    "RdBu",
    "RdYlBu",
    "RdYlGn",
    "Spectral",
    "coolwarm",
    "bwr",
    "seismic",
]

uniform_cmaps = ["viridis", "plasma", "inferno", "magma", "cividis"]

#################################################################
###
###         PLOTTING FUNCTIONS
###
#################################################################


def plot_anomaly(product, variable, year, month, cmap, boundz, lta_period):
    """Do an anomaly (e.g. z-score) plot of a variable in a given
    product.

    Parameters
    ----------
    product : str
        One of the supported products: `TAMSAT`, `ERA5` or `MODIS`.
    variable : str
        A reasonable variable name
    year: int
        The year... D'oh
    month: int
        The month in the current year to use for calculations
    cmap : str
        A matplotlib colormap name
    boundz : iter
        A two element iterable with the min-max values for the colormap scale
    lta_period : str
        What LTA period
    """

    ds = get_one_year(product, variable, year)

    vmin, vmax = boundz
    m, s = get_climatology(product, variable, period=lta_period)
    z_score = calculate_z_score(
        ds, curr_month_number=month, clim_mean=m, clim_std=s
    )
    title = f"Anomaly map {product}-{variable} vs LTA {lta_period}"

    if product in ["ERA5", "TAMSAT"]:
        fig = do_map(
            z_score,
            title,
            r"$z$-score",
            contour=True,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    else:
        fig = do_map(
            z_score,
            title,
            r"$z$-score",
            contour=False,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    for ext in ["png", "pdf"]:
        print(
            f"Saving to {product}_anom_{variable}_{lta_period}_zscore.{ext}"
        )
        fig.savefig(
            f"{product}_anom_{variable}_{lta_period}_zscore.{ext}",
            dpi=175,
            bbox_inches="tight",
        )


def plot_field(product, variable, year, month, cmap, boundz):
    """Do an anomaly (e.g. z-score) plot of a variable in a given
    product.

    Parameters
    ----------
    product : str
        One of the supported products: `TAMSAT`, `ERA5` or `MODIS`.
    variable : str
        A reasonable variable name
    year: int
        The year... D'oh
    month: int
        The month in the current year to use for calculations
    cmap : str
        A matplotlib colormap name
    boundz : iter
        A two element iterable with the min-max values for the colormap scale
    lta_period : str
        What LTA period
    """
    ds = get_one_year(product, variable, year)
    title = f"{product}-{variable}"
    units = PRODUCT_UNITS[product][variable]
    if variable in ["precip", "runoff", "rfe_filled"]:
        # This are aggregated monthly fluxes
        curr_month = (
            ds.resample(time="1MS")
            .sum()
            .sel({"time": f"{year}-{month:02d}-01"})
        )
    else:
        # This are averaged monthly fluxes
        curr_month = (
            ds.resample(time="1MS")
            .mean()
            .sel({"time": f"{year}-{month:02d}-01"})
        )
    vmin, vmax = (
        curr_month.chunk({"x": -1, "y": -1})
        .quantile([boundz[0] / 100, boundz[1] / 100])
        .values
    )
    if product in ["ERA5", "TAMSAT"]:
        fig = do_map(
            curr_month,
            title,
            units,
            contour=True,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    else:
        fig = do_map(
            curr_month,
            title,
            units,
            contour=False,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
        )
    for ext in ["png", "pdf"]:
        print(f"Saving to {product}_{variable}.{ext}")
        fig.savefig(
            f"{product}_{variable}.{ext}",
            dpi=175,
            bbox_inches="tight",
        )


#################################################################
###
###         GUI FUNCTION
###
#################################################################


def do_plots_gui():
    # year = dt.datetime.now().year
    year = 2020
    last_dates = check_dates()
    product = "TAMSAT"
    product_sel = widgets.RadioButtons(
        options=["TAMSAT", "ERA5", "MODIS"],
        value="TAMSAT",
        description="Product family",
    )
    variable_sel = widgets.Dropdown(
        options=variable_lists[product], description="Variable"
    )

    months_this_year = widgets.Select(
        options=range(1, last_dates[product] + 1),
        description="Month current year",
    )
    anomaly_calc = widgets.Checkbox(
        value=True, description="Do anomaly calculations"
    )
    sel_cmap = widgets.Dropdown(options=diverging_cmaps, value="Spectral")
    sel_boundz = widgets.FloatRangeSlider(
        value=[-2.5, 2.5],
        min=-6,
        max=6,
        step=0.1,
        description="Scale for anomaly colormap",
    )
    sel_period = widgets.RadioButtons(
        options=periods, description="LTA temporal window", disabled=False
    )

    def on_value_change(change):
        product = change.new
        variable_sel.options = variable_lists[product]

    def on_value_change2(change):
        product = change.new
        months_this_year.options = range(1, last_dates[product] + 1)

    product_sel.observe(on_value_change, "value")
    product_sel.observe(on_value_change2, "value")

    def on_value_change3(change):
        anomaly = change.new
        if anomaly:
            sel_cmap.options = diverging_cmaps
        else:
            sel_cmap.options = uniform_cmaps

    def on_value_change4(change):
        anomaly = change.new
        if anomaly:
            sel_boundz.min = -6
            sel_boundz.max = 6
            sel_boundz.value = [-2.5, 2.5]
            sel_boundz.step = 0.1
        else:
            sel_boundz.min = 0
            sel_boundz.max = 100
            sel_boundz.value = [5, 95]
            sel_boundz.step = 1

    def on_value_change5(change):
        anomaly = change.new
        if anomaly:
            sel_period.layout.visibility = "visible"
        else:
            sel_period.layout.visibility = "hidden"

    anomaly_calc.observe(on_value_change3, "value")
    anomaly_calc.observe(on_value_change4, "value")
    anomaly_calc.observe(on_value_change5, "value")

    def do_plots(**kwds):
        if kwds["anomaly"]:
            plot_anomaly(
                kwds["product"],
                kwds["variable"],
                year,
                kwds["month"],
                kwds["cmap"],
                kwds["boundz"],
                kwds["lta_period"],
            )
        else:
            plot_field(
                kwds["product"],
                kwds["variable"],
                year,
                kwds["month"],
                kwds["cmap"],
                kwds["boundz"],
            )

    widgets.interact_manual(
        do_plots,
        anomaly=anomaly_calc,
        product=product_sel,
        variable=variable_sel,
        cmap=sel_cmap,
        boundz=sel_boundz,
        month=months_this_year,
        lta_period=sel_period,
    )
