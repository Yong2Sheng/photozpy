#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 18 01:15:36 2026

@author: shengyong
"""

from .preset_filters import SDSS_PRIMED
from .telescope import Telescope
from .preset_ccds import RM_Andor, CT_Andor, CT_FLI

SARA_RM = Telescope(
    ccds=RM_Andor,
    filter_sets=SDSS_PRIMED,
)

SARA_CT = Telescope(
    ccds=(CT_Andor, CT_FLI),
    filter_sets=SDSS_PRIMED,
)