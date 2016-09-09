import os
import datetime
from glob import glob
import numpy as np
import pandas as pd
import engarde.decorators as ed

data_path = os.path.dirname(__file__)

_daily_dtypes = dict(tmin=float, tmax=float, prcp=float)
_daily_ranges = dict(tmin=(-100, 100), tmax=(-100, 100), prcp=(0, 1000))


def load(data_set):
    if data_set.lower() in ['seattle', 'telluride']:
        return _load_daymet_data_set(data_set)
    elif data_set.lower() in ['stehekin_daily']:
        return _load_stehekin_daily()
    elif data_set.lower() in ['stehekin_disagg']:
        return _load_stehekin_disagg()
    else:
        return ValueError('Unknown test dataset %s' % data_set)


def _parse_daymet_dates(year, yday):
    return np.array([datetime.datetime(int(y), 1, 1) +
                     datetime.timedelta(days=int(d) - 1)
                     for y, d in zip(year, yday)])


def _rational_daily(df):
    '''Check that at the daily inputs are rational'''
    r = df.tmax > df.tmin
    return r


@ed.has_dtypes(items=_daily_dtypes)
@ed.verify_all(_rational_daily)
@ed.none_missing()
@ed.within_range(_daily_ranges)
def _read_bin(filename, rec_dtype, mults, index):
    return pd.DataFrame(np.fromfile(filename, dtype=rec_dtype).astype(float),
                        index=index) / mults


@ed.has_dtypes(items=_daily_dtypes)
@ed.verify_all(_rational_daily)
@ed.none_missing()
@ed.within_range(_daily_ranges)
def _load_daymet_data_set(data_set='Seattle'):
    if data_set.lower() == 'seattle':
        fname = os.path.join(data_path, 'data/daymet_seattle.csv')
    elif data_set.lower() == 'telluride':
        fname = os.path.join(data_path, 'data/daymet_telluride.csv')
    renames = {}
    units = {}
    df = pd.read_csv(fname, skiprows=6,
                     index_col='datetime',
                     parse_dates={'datetime': ['year', 'yday']},
                     date_parser=_parse_daymet_dates)
    for c in df.columns:
        n, u = c.split(' ', 1)
        renames[c] = n
        units[n] = u

    df = df.rename(columns=renames)

    return df


def _load_stehekin_daily():
    data = {}
    files = glob('./data/vic42_daily/data_*')

    dt = np.dtype([('prcp', 'uint16'), ('tmax', 'int16'),
                   ('tmin', 'int16'), ('wind', 'int16')])
    mults = np.array([40, 100, 100, 100])
    index = pd.date_range(start='1949-01-01', freq='D', periods=20819)

    for fname in files:
        coord = tuple(map(float, os.path.split(fname)[-1].split('_')[-2:]))
        data[coord] = _read_bin(fname, dt, mults, index=index)

    return data


def _load_stehekin_disagg():
    data = {}
    files = glob('./data/vic42_disagg/full_data_*')

    index = pd.date_range(start='1949-01-01 00:00:00', freq='H', periods=26280)

    for fname in files:
        coord = tuple(map(float, os.path.split(fname)[-1].split('_')[-2:]))
        data[coord] = pd.read_table(fname, skiprows=6, index=index,
                                    names=['prcp', 'air_temp', 'shortwave',
                                           'longwave', 'density', 'pressure',
                                           'vp', 'wind'])

    return data
