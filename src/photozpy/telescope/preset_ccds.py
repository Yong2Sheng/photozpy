#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 01:23:39 2026

@author: shengyong
"""

from datetime import date

import astropy.units as u

from .ccd import CCD

RM_Andor = CCD(
    name="RM Andor Ikon-L",
    gain=1.0 * u.electron / u.adu,
    rdnoise=6.3 * u.electron / u.pix,
    valid_from=date(2016, 1, 1),
    native_pixel_scale=0.34 * u.arcsec / u.pix,
)

CT_Andor = CCD(
    name="CT Andor Ikon-L",
    gain=1.0 * u.electron / u.adu,
    rdnoise=7.1 * u.electron / u.pix,
    valid_from=date(2021, 1, 1),
    native_pixel_scale=0.34 * u.arcsec / u.pix,
)

CT_FLI = CCD(
    name="CT FLI",
    gain=2.0 * u.electron / u.adu,
    rdnoise=(9.7 * 2.0) * u.electron / u.pix,
    valid_from=date(2015, 1, 1),
    valid_to=date(2020, 12, 31),
    native_pixel_scale=0.61 * u.arcsec / u.pix,
)
