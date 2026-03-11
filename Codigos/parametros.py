"""
parametros:
    Se encarga de tomar los archivos y calcular
    diversos parametros como la Temperatura y demas
    
Alguna libreria de parametros util
https://sites.uni.edu/morgans/astro/course/Notes/section2/spectraltemps.html

"""


import re
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.visualization import quantity_support
from pathlib import Path
from scipy.signal import convolve
from scipy.signal import find_peaks, peak_prominences, peak_widths


#%% Temperatura

def get_Temp(lamb, flux): #Saca temperatura con Ley Wien
    """
    Ley de Wien: 
        lamb_max*T= b
    """

    b= 2.8976e7 #Angstrom
    index = np.where(flux == np.max(flux))
    """
    Index toma la posicion del maximo del flujo y se la guarda para luego
    poder sacar la longitud de onda que tenga esa posicion en el array de longitudes
    de onda
    """

    lamb_max = lamb[index]

    T= np.round(b/lamb_max)

    return T

def picos (wave, flux, lineas, dist=10):
    '''
    La salida de la función es del tipo [picos,[prom,ancho]]; siendo cada elemento un array. En caso de no encontrar pico devuelve [0,[0,0]]
    '''
    peaks, props = [], []   #se inicilizan las listas vacías (sé que el append es lento de cojones pero macho que son 4 cosas
    for line in lineas:     
        center = np.searchsorted(wave, lineas[line], side='right') #pilla el i para que el i-ésimo componente de los datos sea el más cercano a la linea dada
        a = center - dist    
        b = center + dist
        rango = -1*flux[a:b] #se delimita el rango de búsqueda de picos 
        ind, prop = find_peaks(rango, distance=2*dist) #el find peaks busca el pico más `significativo` que hay en el rango 
        if ind == []:    #por si no hay picos
            ind = 0
            propiedades = [0,0]
        else:    
            ind += a    #se recentra el índice para los datos de longitud de onda dados
            propiedades = [peak_prominences(rango, ind)[0], peak_widths(rango, ind)[0]] #información adicional de los picos
            peaks.append(wave[ind])
            props.append(propiedades) #se añaden los datos a la lista. se pude hacer sin el apend pero es más coñazo porque la entrada es un diccionario
    return peaks, props
    
      
