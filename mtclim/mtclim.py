'''
mtclim: Mountain Climate Simulator
'''

# Mountain Climate Simulator, meteorological forcing disaggregator
# Copyright (C) 2015  Joe Hamman

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import pandas as pd
import numpy as np
from warnings import warn


constants = \
    {'SECPERRAD': 13750.9871,  # seconds per radian of hour angle
     'RADPERDAY': 0.017214,  # radians of Earth orbit per julian day
     'RADPERDEG': 0.01745329,  # radians per degree
     'MINDECL': -0.4092797,  # minimum declination (radians)
     'DAYSOFF': 11.25,  # julian day offset of winter solstice
     'SRADDT': 30.0,  # timestep for radiation routine (seconds)
                      # Note:  Make sure that 3600 % SRADDT == 0
     'MA': 28.9644e-3,  # (kg mol-1) molecular weight of air
     'R': 8.3143,  # (m3 Pa mol-1 K-1) gas law constant
     'G_STD': 9.80665,  # (m s-2) standard gravitational accel.
     'P_STD': 101325.0,  # (Pa) standard pressure at 0. m elevation
     'T_STD': 288.15,  # (K) standard temp at 0. m elevation
     'CP': 1010.0,  # (J kg-1 K-1) specific heat of air
     'LR_STD': 0.0065,  # (-K m-1) standard temperature lapse rate
     # optical airmass by degrees
     'OPTAM': np.array([2.90, 3.05, 3.21, 3.39, 3.69, 3.82, 4.07, 4.37,
                       4.72, 5.12, 5.60, 6.18, 6.88, 7.77, 8.90, 10.39,
                       12.44, 15.36, 19.79, 26.96, 30.00]),
     'KELVIN': 273.15,
     'EPS': 0.62196351,
     }

default_parameters = \
    {'TDAYCOEF': 0.45,  # (dim) daylight air temperature coefficient

     # parameters for the snowpack algorithm
     'SNOW_TCRIT': -6.0,  # (deg C) critical temperature for snowmelt
     'SNOW_TRATE': 0.042,  # (cm/degC/day) snowmelt rate

     # parameters for the radiation algorithm
     # (dim) stands for dimensionless values
     'TBASE': 0.870,  # (dim) max inst. trans., 0m, nadir, dry atm
     'ABASE': -6.1e-5,  # (1/Pa) vapor pressure effect on
                        # transmittance
     'C': 1.5,  # (dim) radiation parameter
     'B0': 0.031,  # (dim) radiation parameter
     'B1': 0.201,  # (dim) radiation parameter
     'B2': 0.185,  # (dim) radiation parameter
     'RAIN_SCALAR': 0.75,  # (dim) correction to trans. for rain day
     'DIF_ALB': 0.6,  # (dim) diffuse albedo for horizon correction
     'SC_INT': 1.32,  # (MJ/m2/day) snow correction intercept
     'SC_SLOPE': 0.096,  # (MJ/m2/day/cm) snow correction slope
     'site_elev': 0.,
     'base_elev': 0.,
     'tmax_lr': 0.0065,
     'tmin_lr': 0.0065,
     'site_isoh': None,
     'base_isoh': None,
     'site_lat': 0.,
     'site_slope': 0.,
     'site_aspect': 0.,
     'site_east_horiz': 0.,
     'site_west_horiz': 0.
     }

default_options = \
    {'SW_PREC_THRESH': 0.,
     'VP_ITER': 'VP_ITER_ALWAYS',
     'MTCLIM_SWE_CORR': False,
     'LW_CLOUD': 'LW_CLOUD_DEARDORFF',
     'LW_TYPE': 'LW_PRATA'}


class MtClim(object):
    '''The University of Montana Mountain Climate Simulator'''

    def __init__(self, data=None, parameters=None, options=None):
        '''
        Initialize MtClim object.

        Parameters
        ----------
        data : pandas.DataFrame, optional
            Input data: pandas DataFrame with at least `tmax`, `tmin`, and
            `prec`. Timestep freq should be `D`.
        parameters : dict-like, optional
            Dictionary of parameters to use apart from the default parameters.
        options : dict-like, optional
            Dictionary of options to use apart from the default options.
        '''

        # Set model parameters
        self.parameters = default_parameters
        if parameters is not None:
            for p, val in parameters.items():
                if p in self.parameters:
                    self.parameters[p] = val
                else:
                    raise ValueError(
                        'Parameter: %s is not a valid parameter' % p)

        # Set model options
        self.options = default_options
        if options is not None:
            for p, val in options.items():
                if p in self.options:
                    self.options[p] = val
                else:
                    raise ValueError(
                        'Option: %s is not a valid option' % p)

        # Initialize data attribute
        if data is not None:
            self.set_data(data)
        else:
            self.data = None

    def set_data(self, data):
        '''set data attribute

        data : pandas.DataFrame
            Input data: pandas DataFrame with at least `tmax`, `tmin`, and
            `prec`. Timestep freq should be `D`.
        '''
        assert type(data) == pd.DataFrame
        assert all([v in data for v in ['tmax', 'tmin', 'prcp']])
        self.data = data.resample('D', how='mean')
        self.ndays = (data.index[-1] - data.index[0]).days

    def __repr__(self):
        r = 'MtClim object\nparameters: {0}\noptions: {1}\ndata: {2}'.format(
            self.parameters, self.options, self.data.head())
        return r

    def __str__(self):
        return 'MtClim object'

    def calc_tair(self):
        '''
        Calculates daily air temperatures.
        '''
        # calculate elevation difference in meters
        dz = (self.parameters['site_elev'] -
              self.parameters['base_elev'])

        # apply lapse rate corrections to tmax and tmin

        # Since tmax lapse rate usually has a larger absolute value than tmin
        # lapse rate, it is possible at high elevation sites for these
        # corrections to result in tmin > tmax. Check for that occurrence and
        # force tmin = corrected tmax - 0.5 deg C.
        self.data['s_tmax'] = self.data['tmax'] + \
            (dz * self.parameters['tmax_lr'])
        self.data['s_tmin'] = self.data['tmin'] + \
            (dz * self.parameters['tmin_lr'])

        self.data['s_tmin'] = np.minimum(self.data['s_tmin'],
                                         self.data['s_tmax'] - 0.5)

        # derived temperatures
        tmean = self.data[['s_tmax', 's_tmin']].mean(axis=1)
        self.data['s_tday'] = ((self.data['s_tmax'] - tmean) *
                               self.parameters['TDAYCOEF']) + tmean

    def calc_prcp(self):
        '''
        Calculates daily total precipitation

        Tests to write:
        - catch TypeError when isohs are None.
        - no change when isohs are not set
        '''
        self.data['s_prcp'] = self.data['prcp']
        try:
            self.data['s_prcp'] *= (self.parameters['site_isoh'] /
                                    self.parameters['base_isoh'])
        except TypeError:
            if (self.parameters['site_isoh'] is None and
                    self.parameters['base_isoh'] is None):
                pass
            else:
                raise ValueError('Only one of base_isoh and site_isoh is set.')

    def snowpack(self):
        '''
        estimates the accumulation and melt of snow for radiation algorithm
        corrections

        Tests to write:
        - runs to completion
        '''
        # initialize SWE array
        self.data['s_swe'] = 0.

        # first pass
        self._simple_snowpack(0.)

        # use the first pass to set the initial snowpack conditions for the
        # first day of data
        start_yday = self.data.index.dayofyear[0]
        prev_yday = (self.data.index[0] - pd.Timedelta(1, unit='D')).dayofyear
        count = 0
        swe_sum = 0.
        for i in range(self.ndays):
            if (self.data.index.dayofyear[i] == start_yday or
                    self.data.index.dayofyear[i] == prev_yday):
                count += 1
                swe_sum += self.data['s_swe'][i]

        # Proceed with correction if there are valid days to reinitialize
        # the snowpack estiamtes. Otherwise use the first-pass estimate.
        if count:
            snowpack = swe_sum / count
            self._simple_snowpack(snowpack)

    def _simple_snowpack(self, snowpack):
        '''
        Tests to write:
        - runs to completion
        - different results when snowpack starts is 0/nonzero
        '''
        for i in range(self.ndays):
            if (self.data['s_tmin'][i] <= self.parameters['SNOW_TCRIT']):
                snowpack += self.data['s_prcp'][i]
            else:
                snowpack -= self.parameters['SNOW_TRATE'] * \
                    (self.data['s_tmin'][i] - self.parameters['SNOW_TCRIT'])
            snowpack = np.maximum(snowpack, 0.)
            self.data['s_swe'][i] = snowpack

    def calc_srad_humidity_iterative(self, tol=0.01, win_type='boxcar'):
        '''
        Iterative estimation of shortwave radiation and humidity'''
        ndays = self.ndays

        daylength = np.zeros(366)
        window = np.zeros(ndays + 90)
        tinystepspday = 86400 / constants['SRADDT']
        tiny_radfract = np.zeros((366, tinystepspday))
        ttmax0 = np.zeros(366)
        flat_potrad = np.zeros(366)
        slope_potrad = np.zeros(366)
        t_fmax = np.zeros(ndays)
        self.data['s_tfmax'] = 0.

        # calculate diurnal temperature range for transmittance calculations
        self.data['tmax'] = np.maximum(self.data['tmax'], self.data['tmin'])
        dtr = self.data['tmax'] - self.data['tmin']

        # smooth dtr array: After Bristow and Campbell, 1984
        # use 30-day antecedent smoothing window
        sm_dtr = pd.rolling_window(dtr, window=30, freq='D',
                                   win_type=win_type).fillna(method='bfill')
        if self.ndays <= 30:
            warn('Timeseries is shorter than rolling mean window, filling '
                 'missing values with unsmoothed data.')
            sm_dtr.fillna(dtr, inplace=True)

        # calculate the annual total precip
        sum_prcp = self.data['s_prcp'].values.sum()

        ann_prcp = (sum_prcp / self.ndays) * 365.25
        if (ann_prcp == 0.):
            ann_prcp = 1.0

        # Generate the effective annual precip, based on a 3-month
        # moving-window. Requires some special case handling for the
        # beginning of the record and for short records.

        # check if there are at least 90 days in this input file, if not,
        # use a simple total scaled to effective annual precip
        if (ndays < 90):
            sum_prcp = self.data['s_prcp'].values.sum()

            effann_prcp = (sum_prcp / self.ndays) * 365.25

            # if the effective annual precip for this period
            # is less than 8 cm, set the effective annual precip to 8 cm
            # to reflect an arid condition, while avoiding possible
            # division-by-zero errors and very large ratios (PET/Pann)
            effann_prcp = np.maximum(effann_prcp, 8.)

            parray = effann_prcp

        else:
            # Check if the yeardays at beginning and the end of this input file
            # match up. If so, use parts of the three months at the end
            # of the input file to generate effective annual precip for
            # the first 3-months. Otherwise, duplicate the first 90 days
            # of the record.
            start_yday = self.data.index.dayofyear[0]
            end_yday = self.data.index.dayofyear[ndays - 1]
            if (start_yday != 1):
                if end_yday == start_yday - 1:
                    isloop = True
            else:
                if end_yday == 365 or end_yday == 366:
                    isloop = True

            # fill the first 90 days of window
            for i in range(90):
                if (isloop):
                    window[i] = self.data['s_prcp'][ndays - 90 + i]
                else:
                    window[i] = self.data['s_prcp'][i]
            # fill the rest of the window array
            window[90:] = self.data['s_prcp']

            # for each day, calculate the effective annual precip from
            # scaled 90-day total
            for i in range(self.ndays):
                sum_prcp = 0.
                for j in range(90):
                    sum_prcp += window[i + j]
                    sum_prcp = (sum_prcp / 90.) * 365.25

                # if the effective annual precip for this 90-day period
                # is less than 8 cm, set the effective annual precip to 8 cm
                # to reflect an arid condition, while avoiding possible
                # division-by-zero errors and very large ratios (PET/Pann)
                if sum_prcp < 8.:
                    parray[i] = sum_prcp

        # start of the main radiation algorithm

        # before starting the iterative algorithm between humidity and
        # radiation, calculate all the variables that don't depend on
        # humidity so they only get done once.

        # STEP (1) calculate pressure ratio (site/reference) = f(elevation)
        t1 = 1.0 - (constants['LR_STD'] * self.parameters['site_elev']) / \
            constants['T_STD']
        t2 = constants['G_STD'] / (constants['LR_STD'] *
                                   (constants['R'] / constants['MA']))
        pratio = np.power(t1, t2)

        # STEP (2) correct initial transmittance for elevation
        trans1 = np.power(self.parameters['TBASE'], pratio)

        # STEP (3) build 366-day array of ttmax0, potential rad, and daylength
        # precalculate the transcendentals
        lat = self.parameters['site_lat']
        # check for (+/-) 90 degrees latitude, throws off daylength calc
        lat *= constants['RADPERDEG']
        if (lat > np.pi / 2.):
            lat = np.pi / 2.
        if (lat < -np.pi / 2.):
            lat = -np.pi / 2.
        coslat = np.cos(lat)
        sinlat = np.sin(lat)
        cosslp = np.cos(self.parameters['site_slope'] * constants['RADPERDEG'])
        sinslp = np.sin(self.parameters['site_slope'] * constants['RADPERDEG'])
        cosasp = np.cos(self.parameters['site_aspect'] *
                        constants['RADPERDEG'])
        sinasp = np.sin(self.parameters['site_aspect'] *
                        constants['RADPERDEG'])
        # cosine of zenith angle for east and west horizons
        coszeh = np.cos(np.pi / 2. - (self.parameters['site_east_horiz'] *
                                      constants['RADPERDEG']))
        coszwh = np.cos(np.pi / 2. - (self.parameters['site_west_horiz'] *
                                      constants['RADPERDEG']))

        # sub-daily time and angular increment information
        dt = constants['SRADDT']  # set timestep
        dh = dt / constants['SECPERRAD']  # calculate hour-angle step

        # begin loop through yeardays
        for i in range(365):
            # calculate cos and sin of declination
            decl = constants['MINDECL'] * np.cos((i + constants['DAYSOFF']) *
                                                 constants['RADPERDAY'])
            cosdecl = np.cos(decl)
            sindecl = np.sin(decl)

            # do some precalculations for beam-slope geometry (bsg)
            bsg1 = -sinslp * sinasp * cosdecl
            bsg2 = (-cosasp * sinslp * sinlat + cosslp * coslat) * cosdecl
            bsg3 = (cosasp * sinslp * coslat + cosslp * sinlat) * sindecl

            # calculate daylength as a function of lat and decl
            cosegeom = coslat * cosdecl
            sinegeom = sinlat * sindecl
            coshss = -(sinegeom) / cosegeom
            if (coshss < -1.0):
                coshss = -1.0  # 24-hr daylight
            if (coshss > 1.0):
                coshss = 1.0  # 0-hr daylight
            hss = np.cos(coshss)  # hour angle at sunset (radians)
            # daylength (seconds)
            daylength[i] = 2.0 * hss * constants['SECPERRAD']

            if (daylength[i] > 86400):
                daylength[i] = 86400

            # solar constant as a function of yearday (W/m^2)
            sc = 1368.0 + 45.5 * np.sin((2.0 * np.pi * i / 365.25) + 1.7)

            # extraterrestrial radiation perpendicular to beam, total over
            # the timestep (J)
            dir_beam_topa = sc * dt

            sum_trans = 0.
            sum_flat_potrad = 0.
            sum_slope_potrad = 0.

            # begin sub-daily hour-angle loop, from -hss to hss
            for h in np.arange(-hss, hss, dh):
                # precalculate cos and sin of hour angle
                cosh = np.cos(h)
                sinh = np.sin(h)

                # calculate cosine of solar zenith angle
                cza = cosegeom * cosh + sinegeom

                # calculate cosine of beam-slope angle
                cbsa = sinh * bsg1 + cosh * bsg2 + bsg3

                # check if sun is above a flat horizon
                if (cza > 0.):
                    # when sun is above the ideal (flat) horizon, do all the
                    # flat-surface calculations to determine daily total
                    # transmittance, and save flat-surface potential radiation
                    # for later calculations of diffuse radiation

                    # potential radiation for this time period, flat surface,
                    # top of atmosphere
                    dir_flat_topa = dir_beam_topa * cza

                    # determine optical air mass
                    am = 1.0 / (cza + 0.0000001)
                    if (am > 2.9):
                        ami = int((np.cos(cza) / constants['RADPERDEG'])) - 69
                        if (ami < 0):
                            ami = 0
                        if (ami > 20):
                            ami = 20
                        am = constants['OPTAM'][ami]

                    # correct instantaneous transmittance for this optical
                    # air mass
                    trans2 = np.power(trans1, am)

                    # instantaneous transmittance is weighted by potential
                    # radiation for flat surface at top of atmosphere to get
                    # daily total transmittance
                    sum_trans += trans2 * dir_flat_topa

                    # keep track of total potential radiation on a flat
                    # surface for ideal horizons
                    sum_flat_potrad += dir_flat_topa

                    # keep track of whether this time step contributes to
                    # component 1 (direct on slope)
                    if ((h < 0. and cza > coszeh and cbsa > 0.) or
                            (h >= 0. and cza > coszwh and cbsa > 0.)):
                        # sun between east and west horizons, and direct on
                        # slope. this period contributes to component 1
                        sum_slope_potrad += dir_beam_topa * cbsa
                else:
                    dir_flat_topa = -1

                tinystep = (12 * 3600 + h * constants['SECPERRAD']) / \
                    constants['SRADDT']
                if (tinystep < 0):
                    tinystep = 0
                if (tinystep > tinystepspday - 1):
                    tinystep = tinystepspday - 1
                if (dir_flat_topa > 0):
                    tiny_radfract[i][tinystep] = dir_flat_topa
                else:
                    tiny_radfract[i][tinystep] = 0

            if (daylength[i] and sum_flat_potrad > 0):
                tiny_radfract[i] /= sum_flat_potrad

            # calculate maximum daily total transmittance and daylight average
            # flux density for a flat surface and the slope
            if (daylength[i]):
                ttmax0[i] = sum_trans / sum_flat_potrad
                flat_potrad[i] = sum_flat_potrad / daylength[i]
                slope_potrad[i] = sum_slope_potrad / daylength[i]
            else:
                ttmax0[i] = 0.
                flat_potrad[i] = 0.
                slope_potrad[i] = 0.

        # force yearday 366 = yearday 365
        ttmax0[365] = ttmax0[364]
        flat_potrad[365] = flat_potrad[364]
        slope_potrad[365] = slope_potrad[364]
        daylength[365] = daylength[364]

        tiny_radfract[365] = tiny_radfract[364]

        # STEP (4)  calculate the sky proportion for diffuse radiation

        # uses the product of spherical cap defined by average horizon angle
        # and the great-circle truncation of a hemisphere. this factor does not
        # vary by yearday.
        avg_horizon = (self.parameters['site_east_horiz'] +
                       self.parameters['site_west_horiz']) / 2.0
        horizon_scalar = 1.0 - np.sin(avg_horizon * constants['RADPERDEG'])
        if (self.parameters['site_slope'] > avg_horizon):
            slope_excess = self.parameters['site_slope'] - avg_horizon
        else:
            slope_excess = 0.
        if (2.0 * avg_horizon > 180.):
            slope_scalar = 0.
        else:
            slope_scalar = 1.0 - (slope_excess / (180.0 - 2.0 * avg_horizon))
            if (slope_scalar < 0.):
                slope_scalar = 0.

        sky_prop = horizon_scalar * slope_scalar

        # b parameter, and t_fmax not varying with Tdew, so these can be
        # calculated once, outside the iteration between radiation and humidity
        # estimates. Requires storing t_fmax in an array.
        # b parameter from 30-day average of DTR
        b = self.parameters['B0'] + self.parameters['B1'] * \
            np.exp(-self.parameters['B2'] * sm_dtr)

        # proportion of daily maximum transmittance
        t_fmax = 1.0 - 0.9 * np.exp(-b * np.power(dtr, self.parameters['C']))

        # correct for precipitation if this is a rain day
        inds = np.nonzero(self.data['prcp'] >
                          self.options['SW_PREC_THRESH'])[0]
        t_fmax[inds] *= self.parameters['RAIN_SCALAR']
        self.data['s_tfmax'] = t_fmax

        # Initial values of vapor pressure, etc
        if 'tdew' in self.data:
            # Observed Tdew supplied
            tdew = self.data['tdew']
        else:
            # Estimate Tdew
            tdew = self.data['s_tmin']
        if 's_hum' in self.data:
            # Observed vapor pressure supplied
            pva = self.data['s_hum']
        else:
            # convert dewpoint to vapor pressure
            pva = svp(tdew)

        # Other values needed for srad_humidity calculation
        pa = atm_pres(self.parameters['site_elev'])
        yday = self.data.index.dayofyear - 1
        self.data['s_dayl'] = daylength[yday]
        tdew_save = tdew
        pva_save = pva

        # Initial estimates of solar radiation, cloud fraction, etc.
        tdew, pva, pet = self._compute_srad_humidity_onetime(
            tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop, daylength,
            parray, pa, dtr)

        # estimate annual PET
        sum_pet = pet.values.sum()
        ann_pet = (sum_pet / self.ndays) * 365.25

        # Reset humidity terms if no iteration desired
        if (('tdew' in self.data) or ('s_hum' in self.data) or
                (self.options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
                 ann_pet / ann_prcp >= 2.5)):
                tdew = tdew_save[:]
                pva = pva_save[:]

        # Set up srad-humidity iterations
        if (self.options['VP_ITER'].upper() == 'VP_ITER_ALWAYS' or
            (self.options['VP_ITER'].upper() == 'VP_ITER_ANNUAL' and
             ann_pet / ann_prcp >= 2.5) or
                self.options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
            if (self.options['VP_ITER'].upper() == 'VP_ITER_CONVERGE'):
                max_iter = 100
            else:
                max_iter = 2
        else:
            max_iter = 1

        # srad-humidity iterations
        iter_i = 1
        rmse_tdew = tol + 1
        while (rmse_tdew > tol and iter_i < max_iter):
            tdew_save = tdew[:]

            tdew, pva, pet = self._compute_srad_humidity_onetime(
                tdew, pva, ttmax0, flat_potrad, slope_potrad, sky_prop,
                daylength, parray, pa, dtr)
            rmse_tdew = 0
            for i in range(self.ndays):
                # use rmse function and vectorize
                rmse_tdew += (tdew[i] - tdew_save[i]) * \
                    (tdew[i] - tdew_save[i])
                rmse_tdew /= self.ndays
            rmse_tdew = np.power(rmse_tdew, 0.5)
            iter_i += 1

        # save humidity in output data structure
        if 's_hum' not in self.data:
            self.data['s_hum'] = pva
        # output humidity as vapor pressure deficit (Pa)
        # calculate saturated VP at tday
        pvs = svp(self.data['s_tday'])
        vpd = pvs - pva
        self.data['s_vpd'] = np.maximum(vpd, 0.)

    def _compute_srad_humidity_onetime(self, tdew, pva, ttmax0, flat_potrad,
                                       slope_potrad, sky_prop, daylength,
                                       parray, pa, dtr):
        '''
        Initial estimates of solar radiation, cloud fraction, etc.

        Parameters
        ----------
        tdew : pandas.Series
            description
        pva : pandas.Series
            description
        ttmax0 : pandas.Series
            description
        flat_potrad : pandas.Series
            description
        slope_potrad : pandas.Series
            description
        sky_prop : pandas.Series
            description
        daylength : pandas.Series
            description
        parray : pandas.Series
            description
        pa : pandas.Series
            description
        dtr : pandas.Series
            description

        Returns
        ----------
        tdew : pandas.Series
            description
        pva : pandas.Series
            description
        pet : pandas.Series
            description
        '''

        yday = self.data.index.dayofyear - 1

        # Compute SW radiation
        t_tmax = ttmax0[yday] + (self.parameters['ABASE'] * pva)

        # this is mainly for the case of observed VP supplied, for
        # which t_tmax sometimes ends up being negative (when potential
        # radiation is low and VP is high)
        t_tmax = np.minimum(t_tmax, 0.0001)

        self.data['s_ttmax'] = t_tmax

        # final daily total transmittance
        t_final = t_tmax * self.data['s_tfmax']
        print(t_tmax)
        print(self.data['s_tfmax'])
        assert not np.isnan(t_final.values).any()

        # estimate fraction of radiation that is diffuse, on an
        # instantaneous basis, from relationship with daily total
        # transmittance in Jones (Plants and Microclimate, 1992)
        # Fig 2.8, p. 25, and Gates (Biophysical Ecology, 1980)
        # Fig 6.14, p. 122.
        pdif = np.clip(-1.25 * t_final + 1.25, 0., 1.)

        # estimate fraction of radiation that is direct, on an
        # instantaneous basis
        pdir = 1.0 - pdif

        # the daily total radiation is estimated as the sum of the
        # following two components:
        # 1. The direct radiation arriving during the part of
        # the day when there is direct beam on the slope.
        # 2. The diffuse radiation arriving over the entire daylength
        # (when sun is above ideal horizon).

        # component 1
        srad1 = slope_potrad[yday] * t_final * pdir

        # component 2 (diffuse)
        # includes the effect of surface albedo in raising the diffuse
        # radiation for obstructed horizons
        srad2 = (flat_potrad[yday] * t_final * pdif) * \
            (sky_prop + self.parameters['DIF_ALB'] * (1.0 - sky_prop))

        # snow pack influence on radiation
        sc = np.zeros_like(self.data['s_swe'])
        if (self.options['MTCLIM_SWE_CORR']):
            inds = np.nonzero(self.data['s_swe'] > 0. * daylength[yday] > 0.)
            # convert to W/m2 and check for zero daylength
            # snow correction in J/m2/day
            sc[inds] = (1.32 + 0.096 * self.data['s_swe'][inds]) *\
                1.0e6 / daylength[yday][inds]

            # set a maximum correction of 100 W/m2
            sc = np.maximum(sc, 100.)  # JJH - this is fishy
            # this could also be sc[inds] = np.maximum(sc[inds], 100.) optimize

        # save daily radiation
        # save cloud transmittance when rad is an input
        if 's_srad' in self.data:
            potrad = (srad1 + srad2 + sc) * daylength[yday] / t_final / 86400
            self.data['s_tfmax'] = 1.0
            inds = np.nonzero((potrad > 0.) * (self.data['s_srad'] > 0.) *
                              (daylength[yday] > 0))[0]
            # both of these are 24hr mean rad. here
            self.data['s_tfmax'][inds] = (self.data['s_srad'][inds] /
                                          (potrad[inds] * t_tmax[inds]))
            self.data['s_tfmax'] = np.maximum(self.data['s_tfmax'], 1.)
        else:
            self.data['s_srad'] = srad1 + srad2 + sc

        if (self.options['LW_CLOUD'].upper() == 'LW_CLOUD_DEARDORFF'):
            self.data['s_tskc'] = (1. - self.data['s_tfmax'])
        else:
            self.data['s_tskc'] = np.sqrt((1. - self.data['s_tfmax']) / 0.65)
        self.data['s_fdir'] = pdir

        # Compute PET using SW radiation estimate, and update Tdew, pva **
        tmink = self.data['s_tmin'] + constants['KELVIN']
        pet = calc_pet(self.data['s_srad'], self.data['s_tday'], pa,
                       self.data['s_dayl'])

        # calculate ratio (PET/effann_prcp) and correct the dewpoint
        ratio = pet / parray
        self.data['s_ppratio'] = ratio * 365.25
        tdewk = tmink * (-0.127 + 1.121 *
                         (1.003 - 1.444 * ratio + 12.312 *
                          np.power(ratio, 2) - 32.766 * np.power(ratio, 3))
                         + 0.0006 * dtr)
        tdew = tdewk - constants['KELVIN']

        pva = svp(tdew)

        return tdew, pva, pet

    def calc_longwave(self):
        '''This routine estimates long wave radiation based on the fractional
        cloud cover (self.data['s_tskc']), the current air temperature (C), and
        the atmospheric vapor pressure (Pa).
        '''

        # See Bras, R. F. , "Hydrology, an introduction to hydrologic science",
        # Addison-Wesley, 1990, p. 42-45

        # convert to Kelvin
        air_temp = self.data['s_tday'] + constants['KELVIN']
        # convert to mbar
        vp = self.data['s_vp'] / 100

        if (self.options['LW_TYPE'].upper() == 'TVA'):
            # TVA (1972) - see Bras 2.35
            emissivity_clear = 0.740 + 0.0049 * vp
        elif (self.options['LW_TYPE'].upper() == 'ANDERSON'):
            # Anderson (1964)
            emissivity_clear = 0.68 + 0.036 * np.power(vp, 0.5)
        elif (self.options['LW_TYPE'].upper() == 'BRUTSAERT'):
            # Brutsaert (1975)
            x = vp / air_temp
            emissivity_clear = 1.24 * np.power(x, 0.14285714)
        elif (self.options['LW_TYPE'].upper() == 'SATTERLUND'):
            # Satterlund (1979)
            emissivity_clear = 1.08 * \
                (1 - np.exp(-1 * np.power(vp, (air_temp / 2016))))
        elif (self.options['LW_TYPE'].upper() == 'IDSO'):
            # Idso (1981)
            emissivity_clear = 0.7 + 5.95e-5 * vp * np.exp(1500 / air_temp)
        elif (self.options['LW_TYPE'].upper() == 'PRATA'):
            # Prata (1996)
            x = 46.5 * vp / air_temp
            emissivity_clear = 1 - (1 + x) * \
                np.exp(-1 * np.power((1.2 + 3 * x), 0.5))

        if (self.options['LW_CLOUD'].upper() == 'CLOUD_DEARDORFF'):
            # Assume emissivity of clouds is 1.0, and that total emissivity is
            # weighted average of cloud emission plus clear-sky emission,
            # weighted by fraction of sky occupied by each
            # (method of Deardorff, 1978)
            cloudfrac = self.data['s_tskc']
            # Deardorff (1978)
            emissivity = cloudfrac * 1.0 + (1 - cloudfrac) * emissivity_clear
        else:
            # see Bras 2.43
            cloudfactor = 1.0 + (0.17 * self.data['s_tskc'] *
                                 self.data['s_tskc'])
            emissivity = cloudfactor * emissivity_clear

        return emissivity * constants['STEFAN_B'] * np.power(air_temp, 4)


def calc_pet(rad, ta, pa, dayl, dt=0.2):
    '''
    calculates the potential evapotranspiration for aridity corrections in
    `calc_vpd()`, according to Kimball et al., 1997

    Parameters
    ----------
    rad : scalar or numpy.ndarray
        daylight average incident shortwave radiation (W/m2)
    ta : scalar or numpy.ndarray
        daylight average air temperature (deg C)
    pa : scalar or numpy.ndarray
        air pressure (Pa)
    dayl : scalar or numpy.ndarray
        daylength (s)
    dt : scalar, optional
        offset for saturation vapor pressure calculation, default = 0.2

    Returns
    ----------
    pet : scalar or numpy.ndarray
        Potential evapotranspiration (cm/day)
    '''
    # rnet       # (W m-2) absorbed shortwave radiation avail. for ET
    # lhvap      # (J kg-1) latent heat of vaporization of water
    # gamma      # (Pa K-1) psychrometer parameter
    # dt = 0.2   # offset for saturation vapor pressure calculation
    # t1, t2     # (deg C) air temperatures
    # pvs1, pvs2 # (Pa)   saturated vapor pressures
    # pet        # (kg m-2 day-1) potential evapotranspiration
    # s          # (Pa K-1) slope of saturated vapor pressure curve

    # calculate absorbed radiation, assuming albedo = 0.2  and ground
    # heat flux = 10% of absorbed radiation during daylight
    rnet = rad * 0.72

    # calculate latent heat of vaporization as a function of ta
    lhvap = 2.5023e6 - 2430.54 * ta

    # calculate the psychrometer parameter: gamma = (cp pa)/(lhvap epsilon)
    # where:
    # cp       (J/kg K)   specific heat of air
    # epsilon  (unitless) ratio of molecular weights of water and air
    gamma = constants['CP'] * pa / (lhvap * constants['EPS'])

    # estimate the slope of the saturation vapor pressure curve at ta
    # temperature offsets for slope estimate
    t1 = ta + dt
    t2 = ta - dt

    # calculate saturation vapor pressures at t1 and t2, using formula from
    # Abbott, P.F., and R.C. Tabony, 1985. The estimation of humidity
    # parameters. Meteorol. Mag., 114:49-56.
    pvs1 = svp(t1)
    pvs2 = svp(t2)

    # calculate slope of pvs vs. T curve near ta
    s = (pvs1 - pvs2) / (t1 - t2)
    # can this be s = svp_slope(ta)? JJH

    # calculate PET using Priestly-Taylor approximation, with coefficient
    # set at 1.26. Units of result are kg/m^2/day, equivalent to mm water/day
    pet = (1.26 * (s / (s + gamma)) * rnet * dayl) / lhvap

    # return a value in centimeters/day, because this value is used in a ratio
    # to annual total precip, and precip units are centimeters
    return (pet / 10.)


def atm_pres(elev):
    '''atmospheric pressure (Pa) as a function of elevation (m)

    Parameters
    ----------
    elev : scalar or numpy.ndarray
        Elevation (meters)

    Returns
    -------
    pressure : scalar or numpy.ndarray
        Atmospheric pressure at elevation `elev` (Pa)

    References
    ----------
    * Iribane, J.V., and W.L. Godson, 1981. Atmospheric Thermodynamics, 2nd
      Edition. D. Reidel Publishing Company, Dordrecht, The Netherlands.
      (p. 168)
    '''
    t1 = 1.0 - (constants['LR_STD'] * elev) / constants['T_STD']
    t2 = constants['G_STD'] / (constants['LR_STD'] * (constants['R'] /
                                                      constants['MA']))

    return constants['P_STD'] * np.power(t1, t2)


def svp(temp, a=0.61078, b=17.269, c=237.3):
    '''Compute the saturated vapor pressure.

    Parameters
    ----------
    temp : numpy.ndarray
        Temperature (degrees Celsius)

    Returns
    ----------
    pressure : numpy.ndarray
        Saturated vapor pressure at temperature `temp` (Pa)

    References
    ----------
    * Maidment, David R. Handbook of hydrology. McGraw-Hill Inc., 1992.
      Equation 4.2.2.
    '''

    svp = a * np.exp((b * temp) / (c + temp))

    inds = np.nonzero(temp < 0.)[0]
    svp[inds] *= 1.0 + .00972 * temp[inds] + .000042 * np.power(temp[inds], 2)

    return svp * 1000.


def svp_slope(temp, a=0.61078, b=17.269, c=237.3):
    '''Compute the gradient of the saturated vapor pressure as a function of
    temperature.

    Parameters
    ----------
    temp : numpy.ndarray
        Temperature (degrees Celsius)

    Returns
    -------
    gradient : numpy.ndarray
        Gradient of d(svp)/dT.

    References
    ----------
    * Maidment, David R. Handbook of hydrology. McGraw-Hill Inc., 1992.
      Equation 4.2.3.
    '''

    return (b * c) / ((c + temp) * (c + temp)) * svp(temp, a=a, b=b, c=c)
