# -*- coding: utf-8 -*-
"""
Load_Data:
    Se encarga de obtener los espectros necesarios
"""
#%% Librerias
import numpy as np
from astropy.io import fits
#%% Funciones

def Load_Dat(filename, path = "Estrellas_Problema", endMessage = False):
    """
    Load_Dat
    
    Dado el nombre de un fichero en formato dat en el path (opcional),
    carga los datos en dos arrays (lambda, flujo)
    """
    base_name = path + "/" + filename
    try:
        lamb, flux = np.loadtxt(base_name,skiprows=3,unpack=True)
    except:
        raise
        print(f"ERROR: No se ha podido cargar correctamente el fichero en path:\n {base_name}")
        return None, None
    if endMessage: 
        print(f"Se ha leido con exito el fichero {filename}. Se ha devuelto en formato [longitud_onda, flujo]")
        print("-"*40)
    return lamb,flux

def Load_Miles(filename,path = "Catalogo_Miles", endMessage = False):
    """
    Load_Miles
    
    Dado un fichero en formato fits, carga los datos en dos arrays (lambda, flujo)
    """
    base_name = path + "/" + filename
    try:
        with fits.open(base_name) as spec1d:
            hdu = spec1d[0]              # o hdul[1] si el espectro está en la primera extensión
            flux = np.ravel(hdu.data)            # array de flujo 1D
            header = hdu.header        # cabecera con las palabras clave WCS
        # Obtengamos las longitudes de onda
        crval1 = header["CRVAL1"]      # valor de longitud de onda en el píxel de referencia[web:4][web:9]
        cdelt1 = header["CDELT1"]      # incremento de longitud de onda por píxel[web:4][web:9]
        crpix1 = header.get("CRPIX1", 1.0)  # píxel de referencia (1\u2011basado); usar 1 si no existe[web:4][web:9].
        # Número de píxeles del espectro
        n_pix = flux.size
        pix = np.arange(n_pix)
        wavelength = ( (pix + 1.0) - crpix1 ) * cdelt1 + crval1
    except:
        print(f"ERROR: NO SE HA PODIDO OBTENER EL ARCHIVO EN EL PATH {base_name}")
        return None, None
    if endMessage:
        print(f"Se ha leido con exito el fichero {filename}. Se ha devuelto en formato [longitud_onda, flujo]")
        print("-"*40)
    return np.array(wavelength,dtype = float), np.array(flux, dtype = float) # Para evitar dtype = float 
