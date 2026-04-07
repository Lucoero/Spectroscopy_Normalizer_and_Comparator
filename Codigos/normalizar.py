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

def Normalizar(lamb,flujo, params = [97],start = 4000, end = -1,  filtro = 'med', iteraciones = 1): #Donde filtro es el tipo de filtro segun filtrar
    """
    Los tipos de filtro posibles son:
        
        med
        mmed
        sg
        h
    El corte es a partir de que longitud de onda (Armstrongs) queremos empezar a normalizar
    ajuste hecho
    flujo normalizado
"""
    ajuste = np.copy(flujo)
    start_index = np.where(lamb > start)[0][0]
    
    normalizada_vieja = flujo[start_index:end]
    lamb_cortado = lamb[start_index:end]
    
    for i in range(iteraciones):
        normalizada_nueva = Filtrar(normalizada_vieja,params, tipo = filtro)
        normalizada_vieja = normalizada_nueva
        
    # Fin ajuste
    ajuste[start_index:end] = normalizada_vieja
    # Nos quitamos el ajuste de 0 porque si no se hacen nans
    ajuste = np.where(ajuste > 0., ajuste, 1.)
    # Ahora normalizamos
    flujo_normalizado = flujo/ajuste
    
    return (lamb_cortado,normalizada_vieja),flujo_normalizado 
    
    
