#!/usr/bin/env jupyter
import numpy as np
import scipy as scp
from scipy.signal import medfilt
from scipy.signal import savgol_filter
from scipy.signal import hilbert


def Filtrar(flux, params, tipo='med'):
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

def Normalizar(flujo, params,corte = 9000 ,  filtro = 'med', iteraciones = 1): #Donde filtro es el tipo de filtro segun filtrar
    """
    Los tipos de filtro posibles son:
        
        med
        mmed
        sg
        h
    El corte es a partir de que longitud de onda (Armstrongs) queremos empezar a normalizar
    ----------
    filtro : TYPE
        DESCRIPTION.
    flujo : TYPE
        DESCRIPTION.
    parametro : TYPE
        DESCRIPTION.
    iteraciones : TYPE
        DESCRIPTION.

    Returns
    -------
    ajuste hecho
    flujo normalizado
"""
#TODO: HACER EL CORTE
    normalizada_vieja = flujo
    for i in range(iteraciones):
        normalizada_nueva = Filtrar(normalizada_vieja,params, tipo= filtro)
        normalizada_vieja = normalizada_nueva
        
    # Ya esta ajustada
    ajuste = normalizada_vieja
    
    # Ahora normalizamos
    flujo_normalizado = flujo/ajuste
    return ajuste,flujo_normalizado
    
    
