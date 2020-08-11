#!/usr/bin/env python

import ipywidgets.widgets as widgets

from .basic_calcs import (
    ERA5_VARIABLES,
    MODIS_VARIABLES,
    TAMSAT_VARIABLES,
    calculate_z_score,
    get_climatology,
    get_era5_ds,
    get_modis_ds,
    get_tamsat_ds,
    plot_z_score,
)

periods = {"2002-2019": "long", "2015-2019": "recent", "2002-2010": "past"}


def plot_tamsat_anomaly():
    widgets.interact_manual(
        tamsat_anomaly,
        variable=widgets.Dropdown(
            options=TAMSAT_VARIABLES, description="TAMSAT variable to plot"
        ),
        lta_period=widgets.Dropdown(
            options=periods, description="LTA temporal window"
        ),
    )


def tamsat_anomaly(variable, lta_period):
    ds = get_tamsat_ds(variable)
    m, s = get_climatology("TAMSAT", variable, period=lta_period)
    z_score = calculate_z_score(ds, clim_mean=m, clim_std=s)
    fig = plot_z_score(z_score, contour=True)
    print(
        f"Saving to tamsat_anom_{variable}_{periods[lta_period]}_zscore.pdf"
    )
    fig.savefig(
        f"tamsat_anom_{variable}_{periods[lta_period]}_zscore.pdf",
        dpi=175,
        bbox_inches="tight",
    )
    fig.savefig(
        f"tamsat_anom_{variable}_{periods[lta_period]}_zscore.png",
        dpi=175,
        bbox_inches="tight",
    )
    print(
        f"Saving to tamsat_anom_{variable}_{periods[lta_period]}_zscore.png"
    )


def plot_era_anomaly():
    widgets.interact_manual(
        era_anomaly,
        variable=widgets.Dropdown(
            options=ERA5_VARIABLES, description="ERA5 variable to plot"
        ),
        lta_period=widgets.Dropdown(
            options=periods, description="LTA temporal window"
        ),
    )


def era_anomaly(variable, lta_period):
    ds = get_era5_ds(variable)
    m, s = get_climatology("ERA", variable, period=lta_period)
    z_score = calculate_z_score(ds, clim_mean=m, clim_std=s)
    fig = plot_z_score(z_score, contour=True)
    print(f"Saving to era5_anom_{variable}_{periods[lta_period]}_zscore.pdf")
    fig.savefig(
        f"era5_anom_{variable}_{periods[lta_period]}_zscore.pdf",
        dpi=175,
        bbox_inches="tight",
    )
    fig.savefig(
        f"era5_anom_{variable}_{periods[lta_period]}_zscore.png",
        dpi=175,
        bbox_inches="tight",
    )
    print(f"Saving to era5_anom_{variable}_{periods[lta_period]}_zscore.png")


def plot_modis_anomaly():
    widgets.interact_manual(
        modis_anomaly,
        variable=widgets.Dropdown(
            options=MODIS_VARIABLES, description="MODIS variable to plot"
        ),
        lta_period=widgets.Dropdown(
            options=periods, description="LTA temporal window"
        ),
    )


def modis_anomaly(variable, lta_period):
    ds = get_modis_ds(product=variable)
    m, s = get_climatology("MODIS", variable, period=lta_period)
    z_score = calculate_z_score(ds, clim_mean=m, clim_std=s)
    fig = plot_z_score(z_score, contour=False)
    print(f"Saving to modis_anom_{variable}_{periods[lta_period]}_zscore.pdf")
    fig.savefig(
        f"modis_anom_{variable}_{periods[lta_period]}_zscore.pdf",
        dpi=175,
        bbox_inches="tight",
    )
    fig.savefig(
        f"modis_anom_{variable}_{periods[lta_period]}_zscore.png",
        dpi=175,
        bbox_inches="tight",
    )
    print(f"Saving to modis_anom_{variable}_{periods[lta_period]}_zscore.png")
