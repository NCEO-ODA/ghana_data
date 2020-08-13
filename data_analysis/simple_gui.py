#!/usr/bin/env python
"""A simple Jupyter notebook GUI to do spatial plots over Ghana.
"""
import datetime as dt

import ipywidgets.widgets as widgets

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
    if product in ["ERA5", "TAMSAT"]:
        fig = do_map(z_score, contour=True, cmap=cmap, vmin=vmin, vmax=vmax)
    else:
        fig = do_map(z_score, contour=False, cmap=cmap, vmin=vmin, vmax=vmax)
    for ext in ["png", "pdf"]:
        print(f"Saving to {product}_anom_{variable}_{lta_period}_zscore{ext}")
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
    if product in ["precip", "runoff", "rfe_filled"]:
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
    vmin, vmax = curr_month.quantile(
        [boundz[0] / 100, boundz[1] / 100]
    ).values
    if product in ["ERA5", "TAMSAT"]:
        fig = do_map(
            curr_month, contour=True, cmap=cmap, vmin=vmin, vmax=vmax
        )
    else:
        fig = do_map(
            curr_month, contour=False, cmap=cmap, vmin=vmin, vmax=vmax
        )
    for ext in ["png", "pdf"]:
        print(f"Saving to {product}_{variable}{ext}")
        fig.savefig(
            f"{product}_{variable}.{ext}", dpi=175, bbox_inches="tight",
        )


#################################################################
###
###         GUI FUNCTION
###
#################################################################


def do_plots_gui():
    year = dt.datetime.now().year
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
