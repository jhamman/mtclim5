import pytest
from . import load_test_data


@pytest.fixture()
def daymet_seattle(scope='module'):
    return load_test_data.load('seattle')


@pytest.fixture()
def daymet_telluride(scope='module'):
    return load_test_data.load('telluride')


@pytest.fixture()
def stehekin_daily(scope='module'):
    return load_test_data.load('stehekin_daily')


@pytest.fixture()
def stehekin_disagg(scope='module'):
    return load_test_data.load('stehekin_disagg')
