MT-CLIM
======
**Mountain Climate Simulator**

## Summary

MT-CLIM is a computer program that uses observations of daily maximum temperature, minimum temperature, and precipitation from one location (the "base") to estimate the temperature, precipitation, radiation, and humidity at another location (the "site"). The base and the site can be at different elevations, and can have different slopes and aspects. Better results are obtained when the base and site are relatively close to one another (at the bottom and the top of a valley, for example).

Temperature estimates at the site are based on the base temperatures and a user-supplied temperature lapse rate. Separate lapse rates can be supplied for daily maximum and minimum temperature. Precipitation estimates at the site are based on the daily record of preciupitation from the base, and a user-specified ratio of annual total precipitation between the site and the base.

The estimation of radiation and humidity are more complex, since these parameters are not assumed to be measured at the base station. Humidity estimates are based on the observation that daily minimum temperature is usually very close to dewpoint temperature. The MT-CLIM algorithm includes corrections to this assumption for arid climates. Radiation estimates are based on the observation that the diurnal temperature range (from minimum temperature to maximum temperature) is closely related to the daily average atmospheric transmittance. In conjunction with information about the latitude, elevation, slope, and aspect of the site, this relationship can be used to estimate daily total radiation with a typical error range of +/- 15%. Some of the limitations of MT-CLIM are the use of a single base station for observations and the need for the user to specify the temperature and precipitation relationships with elevation.

## Development History

The development of MT-CLIM was originally developed by Peter Thornton at the Numerical Terradynamic Simulation Group (NTSG) at the University of Montana.  MT-CLIM has been used as a run-time meteorological forcing disaggregator inside the [Variable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC).

## Version 5

Development of MT-CLIM version 5.0 began with the version of MT-CLIM version 4.3 embedded in the [Variable Infiltration Capacity (VIC) model](https://github.com/UW-Hydro/VIC) version 4.2.  From that point, the C source code was extracted from VIC and converted to [Python](https://www.python.org/).  Version 5.0 will include the following new features:

- [Python](https://www.python.org/)/[Pandas](http://pandas.pydata.org/) interface.
- Algorithmic flexibility (resampling methods, model choices, and parameter space).
- Subaily estimation or temporal disaggregation of all variables.
- Improved solar geometry functions.
- Improved estimation of tranmissivity by calculating the mean cos(solar zenith angle) and integrating that quantity when disaggregating shortwave to sub-daily time steps and calculating transmissivity.
- ASCII text and netCDF input/output options.


