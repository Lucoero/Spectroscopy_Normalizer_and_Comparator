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
    
def categorizar(wave,flux,lineas,cat):
    """
    INDICATIVOS PARA CADA TIPO ESPECTRAL:
        1. TIPO O:
                - Líneas de Halpha,beta y demás algo débiles (por saha se ha ionizado ya todo el H)
                - Importante la linea de He I y He II
            1.1 Para clasificar la luminosidad:
                - Es muy util comparar la profundidad de las lineas He I y He II:
                  A menor T, mayor profundidad de He I (por Saha) y menor de He II (por menor T)
    
        2. TIPO B:
                Aquí hay 3 criterios espectrales importantes:
                1- Se presentan líneas de HeI pero no de HeII, aunque no es tan efectivo con la digitalización de los datos
                2- He I tiene su maximo en B2, y luego van bajando
                3- Este criterio es el preferido SiIV 4089/SiIII 4552 para tipos tempranos, SiIII 4552/SiII 4128
                5- Tambien importa el ratio He I 4471 y Mg II 4481: A menor T, menor He I y mayor Mg II
                6- Para líneas débiles de Si, MgI 4481/HeI 4471 el ratio incrementa con el subtipo hasta que ronda la unidad en B8
        3. TIPO A:
                - Las lineas de H I desaparecen (ojo, estas son las de ionizacion no las de Balmer)
                - Las lineas de He I tambien desaparecen
                - Empiezan a verse lineas de metalicidad. En especial las de Fe II 4383.
                - Sigue habiendo lineas de H8 y H0, y en este tipo acaban alcanzando su maximo.
        4. TIPO F:
                - Se van debilitando las lineas H8, H9
                - Es importante la linea Ca II, que se vuelve super profunda hasta saturar
                - Empieza a haber una fusion entre CaII H + Hepsi
            
        5. TIPO G:
                - Aparece la banda G (molecular). Cerca de los 4320 A
                - La clasificador de luminosidad se basa en comparar la banda G con la de H delta
        6. TIPO K:
                - Importa Ca I 4227 en especial en K5
                - Las bandas moleculares empiezan a dominar
                - Aparece la franja TiO alpha, que es una estructura de balcones entre 4800 y 5600 A
        7. TIPO M:
                - Ca I 4226 es un valle cada vez mas prominente
    """
    return "Hola"
