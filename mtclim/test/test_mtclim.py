import pytest
import pandas as pd
import numpy as np
from io import StringIO

from mtclim.mtclim import MtClim, atm_pres, svp, svp_slope, calc_pet


@pytest.fixture()
def short_daily_dataframe(scope='function'):
    data = u'''tmax tmin prcp
1915-01-01 31.650000 22.420000 1.32
1915-01-02 31.709999 22.459999 5.24
1915-01-03 32.080002 22.900000 0.00
1915-01-04 32.700001 23.049999 4.11
1915-01-05 32.840000 23.030001 0.00'''
    return pd.read_table(StringIO(data), sep=' ', parse_dates=True)


def test_make_instance():
    m = MtClim()
    print(str(m))


@pytest.mark.usefixtures("short_daily_dataframe")
def test_make_instance_with_data(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe)
    print(str(m))


@pytest.mark.usefixtures("short_daily_dataframe")
def test_set_data_method(short_daily_dataframe):
    m = MtClim()
    m.set_data(short_daily_dataframe)


def test_raise_error_with_bad_options():
    with pytest.raises(ValueError):
        MtClim(options={'junk': 0.})


def test_raise_error_with_bad_parameters():
    with pytest.raises(ValueError):
        MtClim(parameters={'junk': 0.})


def test_assert_data_is_dataframe():
    with pytest.raises(AssertionError):
        MtClim(data=np.random.random((19, 3)))


@pytest.mark.usefixtures("short_daily_dataframe")
def test_assert_correct_input_vars(short_daily_dataframe):
    data = short_daily_dataframe.rename(columns={'prcp': 'Precipitation'})
    with pytest.raises(AssertionError):
        MtClim(data=data)


@pytest.mark.usefixtures("short_daily_dataframe")
def test_calc_tair_runs(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe)
    m.calc_tair()


@pytest.mark.usefixtures("short_daily_dataframe")
def test_calc_tair_keeps_tmax_gt_tmin(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe,
               parameters={'tmax_lr': -65, 'site_elev': 3000})
    m.calc_tair()
    assert (m.data['s_tmax'] >= m.data['s_tmin']).all()
    assert (m.data['s_tmax'] >= m.data['s_tday']).all()
    assert (m.data['s_tday'] >= m.data['s_tmin']).all()


@pytest.mark.usefixtures("short_daily_dataframe")
def test_calc_prcp_catch_single_isoh(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe,
               parameters={'site_isoh': None, 'base_isoh': 43})
    with pytest.raises(ValueError):
        m.calc_prcp()


@pytest.mark.usefixtures("short_daily_dataframe")
def test_calc_prcp_isohs_none(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe,
               parameters={'site_isoh': None, 'base_isoh': None})
    m.calc_prcp()


@pytest.mark.usefixtures("short_daily_dataframe")
def test_snowpack(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe)
    m.calc_tair()
    m.calc_prcp()
    m.snowpack()
    assert (m.data['s_swe'] >= 0.).all()


@pytest.mark.usefixtures("short_daily_dataframe")
def test_simple_snowpack(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe)
    m.calc_tair()
    m.calc_prcp()
    m.data['s_swe'] = 0.
    m._simple_snowpack(0.)
    assert (m.data['s_swe'] >= 0.).all()
    m.data['s_swe'] = 0.
    m._simple_snowpack(100)
    assert (m.data['s_swe'] >= 0.).all()


def test_atm_pres():
    p = atm_pres(np.arange(-250, 4000, 50))
    assert type(p) == np.ndarray
    assert (np.diff(p) < 0.).all()


def test_svp():
    temp = np.arange(-20, 20, dtype=np.float)
    p = svp(temp)
    assert type(p) == np.ndarray
    assert (p > 0.).all()


def test_svp_slope():
    temp = np.arange(-20, 20, dtype=np.float)
    ps = svp_slope(temp)
    assert type(ps) == np.ndarray
    assert (ps > 0.).all()


def test_calc_pet():
    nsteps = 50
    rad = np.linspace(0, 400, nsteps)
    ta = np.linspace(-20, 20, nsteps)
    pa = np.linspace(101000, 102000, nsteps)
    dayl = np.linspace(0, 86400, nsteps)
    pet = calc_pet(rad, ta, pa, dayl)
    assert (pet >= 0.).all()
    assert pet[0] == 0.


@pytest.mark.usefixtures("short_daily_dataframe")
def test_calc_srad_humidity_iterative(short_daily_dataframe):
    m = MtClim(data=short_daily_dataframe)
    m.calc_tair()
    m.calc_prcp()
    m.snowpack()
    assert (m.data['s_swe'] >= 0.).all()
    m.calc_srad_humidity_iterative()
