# Extreme Rainfall Occurrence in Brazil — LGCP/SPDE-INLA

Code and processing pipeline for a Bayesian spatio-temporal analysis of extreme rainfall occurrence over Brazil. Daily CHIRPS precipitation (1981–2026) is used to construct spatio-temporal point patterns of threshold exceedances at 100 mm/day and 200 mm/day. Occurrence intensity is modelled as a Log-Gaussian Cox Process with latent trend, seasonal and autoregressive temporal components and a dynamic spatial field represented through the SPDE approach, with inference performed using INLA.

Two specifications are included: a common-component model, in which the temporal structures are shared across the country, and a regionalized model, in which trend, seasonal and autoregressive components are estimated separately for the five Brazilian macroregions (North, Northeast, Center-West, Southeast, South).

Requirements

R (≥ 4.3) with the following packages:

Purpose	Packages
Inference	INLA, fmesher
Spatial handling	sf, sp, terra, deldir, rnaturalearth
Data reading	ncdf4, dplyr, zoo
Figures	ggplot2, patchwork, viridis, ggspatial

INLA is not on CRAN and must be installed from the project repository:

r
install.packages("INLA",
  repos = c(getOption("repos"),
            INLA = "https://inla.r-inla-download.org/R/stable"),
  dep = TRUE)

Data

CHIRPS v2.0 daily precipitation at 0.05° resolution, in NetCDF format, from the Climate Hazards Center:

https://data.chc.ucsb.edu/products/CHIRPS-2.0/global_daily/netcdf/p05

baixachirps.sh downloads one file per year. The raw NetCDF files are not included in this repository — they total several hundred GB and are freely available from the source above.

Pipeline

Scripts run in the order below. Each threshold has its own copy of every script; the 100 mm versions are listed here, and the 200 mm versions are identical apart from the threshold value and file names.

1. baixachirps.sh — downloads the annual CHIRPS NetCDF files.

2. extraipontosr1_100.R — reads each NetCDF file, converts each daily field to a raster, crops and masks it to the Brazilian boundary, and extracts the coordinates and values of all pixels exceeding the threshold. Output: extreme_rain_points100.csv, one row per exceedance pixel per day.

3. createdataset100.R — builds the spatial and temporal structure. Constructs a constrained Delaunay mesh over Brazil (max.edge = c(3, 10), cutoff = 0.8), a one-dimensional temporal mesh over the monthly index, and a Voronoi tessellation of the mesh nodes. Exceedance points are aggregated into counts per Voronoi cell per month, and the exposure term is computed as the area of each cell intersected with the national boundary, multiplied by the temporal weight. Output: modeloextremerain100.Rdata.

4. modelomultitrend100.R — fits the LGCP. Poisson likelihood with the exposure offset, an RW2 long-term trend, an intrinsic seasonal component of period 12, an AR(2) short-memory component, and a Matérn SPDE field with an AR(1) group structure over time. Output: modeloextremerain100intercept.Rdata (common specification) or modeloextremerain100interceptMulti.Rdata (regionalized).

5. grafico100rev.R / graficos100multirev.R — produce the monthly count series, the temporal component decompositions, and the spatial field panels.

Notes on the model specification

The AR(2) component is parameterized through partial autocorrelations (pacf1, pacf2), which is INLA's internal parameterization and guarantees stationarity. The values reported in the parameter tables are therefore partial autocorrelations; the corresponding autoregressive coefficients follow from the Levinson–Durbin relations, φ₁ = ψ₁(1 − ψ₂) and φ₂ = ψ₂, or equivalently via inla.ar.pacf2phi().

The RW2 precision uses a PC prior with param = c(10, 0.1). Remaining hyperparameters use INLA defaults, with initial values set from preliminary runs to speed convergence. Inference uses the empirical Bayes integration strategy (int.strategy = "eb").

Hardware

[Fill in: machine, cores, memory, whether a cluster was used, approximate runtime per model.]

Citation

[Paper citation once accepted.]

License

MIT.
