#!/usr/bin/env python
# -*- coding: utf-8 -*-

__version__ = "0.1.0"

try:
    __TRANSIT_SETUP__
except NameError:
    __TRANSIT_SETUP__ = False

if not __TRANSIT_SETUP__:
    __all__ = ["System", "Central", "Body"]
    from .transit import System, Central, Body
