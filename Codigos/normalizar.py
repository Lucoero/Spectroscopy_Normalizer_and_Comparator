#!/usr/bin/env jupyter
import numpy as np
import scipy as scp
from scipy.signal import medfilt
from scipy.signal import savgol_filter
from scipy.signal import hilbert


def filtro(flux, params, tipo='med'):
    if tipo == 'med':
        return medfilt(flux, *params)
    if tipo == 'mmed':
        s = medfilt(flux, params[0])
        for i in params:
            s = medfilt(s, i)
        return s
    elif tipo == 'sg':
        return savgol_filter(flux, *params)
    elif tipo == 'h':
        return hilbert(flux)
