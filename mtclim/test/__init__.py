import pytest
import os
import importlib

_travis = os.getenv('TRAVIS', False)

_libs = {}
for mod in ['netCDF4', 'xray', 'matplotlib']:
    try:
        importlib.import_module(mod)
        _libs[mod] = True
    except ImportError:
        _libs[mod] = False

requires_netCDF4 = pytest.mark.skipif(_libs['netCDF4'], reason="requires netCDF4")
requires_matplotlib = pytest.mark.skipif(_libs['matplotlib'],
                                         reason="requires matplotlib")
requires_xray = pytest.mark.skipif(_libs['xray'], reason="requires xray")
skip_on_travis = pytest.mark.skipif(_travis, reason="skip on Travis")
