# -*- coding: utf-8 -*-
"""
Load_Data:
    Se encarga de obtener los espectros necesarios
"""
#%% Librerias
import numpy as np
from astropy.io import fits
import tqdm as tqdm
import os as os
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

def Load_Miles(filename,path = "Catalogo_Miles", endMessage = False, returnName = False):
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
        # Y el nombre del bicho
        name = header["OBJECT"]
        # Número de píxeles del espectro
        n_pix = flux.size
        pix = np.arange(n_pix)
        wavelength = ( (pix + 1.0) - crpix1 ) * cdelt1 + crval1
    except:
        print(f"ERROR: NO SE HA PODIDO OBTENER EL ARCHIVO EN EL PATH {base_name}")
        if returnName:
            return None,None,""
        return None, None
    if endMessage:
        print(f"Se ha leido con exito el fichero {filename}. Se ha devuelto en formato [longitud_onda, flujo]")
        print("-"*40)
    if returnName:
        return np.array(wavelength,dtype = float), np.array(flux, dtype = float),name # Para evitar dtype = float 
    return np.array(wavelength,dtype = float), np.array(flux, dtype = float) # Para evitar dtype = float 

def Write_Data(LambFluxLists,NameList,foldPath, showEndMessage = False):
    """
    Write_Data:
        Dada una lista de tuplas de (lambda, Flux), escribe los datos con nombre
        NameList asociado en la carpeta indicada. Usamos formato .dat, no hace falta incluirlo
        en NameList
    """
    n = len(LambFluxLists)
    if not os.path.exists(foldPath):
        print(foldPath + "didnt exist. We created it.")
        os.mkdir(foldPath)
    for i in tqdm.tqdm(range(n)):
        Lamb,Flux = LambFluxLists[i]
        currName = NameList[i]
        try:
            with open(foldPath + "/" + currName + ".dat","w") as file:
                for j in range(len(Lamb)): 
                    file.write(str(Lamb[j]) + "\t\t\t" + str(Flux[j]))
                    file.write("\n")
        except:
            print("ERROR EN LA ESCRITURA DEL FICHERO.")
            pass
    if showEndMessage:
        print(f"Files stored in \n {foldPath} \n")
    return