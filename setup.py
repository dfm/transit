#!/usr/bin/env python
# encoding: utf-8

import sys
import os

try:
    from setuptools import setup, Extension
    setup, Extension
except ImportError:
    from distutils.core import setup, Extension
    setup, Extension

import numpy


if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()


desc = open("README.rst").read()
required = ["numpy"]

include_dirs = [
    "transit",
    numpy.get_include(),
]
src = [
    "transit/_transit.c",
    "transit/quad.cpp",
    "transit/driver.cpp",
]
ext = Extension("transit._transit", src, include_dirs=include_dirs)

setup(
    name="transit",
    version="0.0.1",
    author="Daniel Foreman-Mackey",
    author_email="danfm@nyu.edu",
    packages=["transit"],
    url="http://github.com/dfm/transit",
    license="MIT",
    description="Compute some transits",
    long_description=desc,
    install_requires=required,
    ext_modules=[ext],
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
