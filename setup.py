#!/usr/bin/env python

try:
    from setuptools import setup
except:
    from distutils.core import setup

setup(name='mtclim',
      version='5.0.0',
      description='The University of Montana Mountain Climate Simulator',
      author='Joe Hamman',
      author_email='jhamman@hydro.washington.edu',
      packages=['mtclim'],
      tests_require=['pytest', 'engarde'])
