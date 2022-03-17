
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Replication: The effects of non-pharmaceutical interventions on COVID-19 mortality

This repository provides all materials for replication of results in

> Mader S and Rüttenauer T (2022) The Effects of Non-pharmaceutical
> Interventions on COVID-19 Mortality: A Generalized Synthetic Control
> Approach Across 169 Countries. *Front. Public Health* 10:820642. doi:
> <https://doi.org/10.3389/fpubh.2022.820642>

Date: 2022-03-17

## Set up

The code for replication of the results requires the following folders:
“01\_Script”, “02\_Data”, “03\_Output”. All R Scripts are required in
folder “01\_Script”, all data are required in folder “02\_Data”.

To reproduce the results of the paper, the scripts need to be executed
in order.

The following packages are necessary for reproduction of main results:

``` r
install.packages(colorspace)
install.packages(cowplot)
install.packages(doParallel)
install.packages(ecmwfr)
install.packages(extrafont)
install.packages(feisr)
install.packages(ggplot2)
install.packages(ggridges)
install.packages(grid)
install.packages(gridExtra)
install.packages(gsynth)
install.packages(lfe)
install.packages(mapview)
install.packages(ncdf4)
install.packages(panelView)
install.packages(parallel)
install.packages(plm)
install.packages(scales)
install.packages(splm)
install.packages(texreg)
install.packages(ungeviz)
install.packages(WDI)
```

### Scripts:

-   *01\_Data\_Covid*: Downloads and prepares data from [Our World in
    Data](https://github.com/owid/covid-19-data/tree/master/public/data)
    and
    [OxCGRT](https://www.bsg.ox.ac.uk/research/research-projects/coronavirus-government-response-tracker#data).

-   *02\_Data\_Weather*: Downloads and prepares weather data from [ERA5
    reanalysis](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview).
    **Note: This requires a valid API user and key** (free
    registration)! The downloaded data is large (NetCDF of around
    700MB).

-   *03\_Descriptives*: Produces Descriptive Figures (Online
    Supplement).

-   *04\_Analysis\_daily*: Produces the main results of Figure 1.

-   *05\_Analysis\_daily\_vaccine*: Produces the results on vaccinations
    of Figure 2.

## Data:

All data for this paper are freely available and accessible online. The
sources are documented in the code, and the data is directly accessed in
the code. Nevertheless, for backwards compatibility, we provide the data
used in the paper in the folder “02\_Data”.

## System and version information

Platform: x86\_64-w64-mingw32/x64 (64-bit) and x86\_64-pc-linux-gnu
(HPC)

Version: R version 4.1.0 and 4.0.5 (HPC)
