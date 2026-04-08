#!/usr/bin/env jupyter
import numpy as np
import scipy as scp
from scipy.signal import medfilt
from scipy.signal import savgol_filter
from scipy.signal import hilbert
import os as os
import tqdm as tqdm

import Load_Data as LD

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
    
def Normalise_Folder(objPathFolder,endPathFolder):
    """
    Normalise_Folder:
        Dado un path con espectros en formato fits o dat,
        y un path de salida, cargamos los datos, los normalizamos,
        y los guardamos en ficheros dat con el nombre original
        y el de la estrella asociado
        Notese que respecto a los fits podemos perder otros datos como la 
        luminosidad, o el corrimiento al rojo, etc.
    """
    if not os.path.exists(objPathFolder):
        print(f"ERROR: WE CANT FIND THE SPECTRAS PATH \n{objPathFolder}")
        return None
    FilesArr = os.listdir(objPathFolder)
    n = len(FilesArr)
    print("\n STARTING NORMALIZATION PROCEDURE...")
    LambFluxLists = []
    outNameList = []
    for i in tqdm.tqdm(range(n)):
        currFile = FilesArr[i]
        lamb, flux,starName = LD.Load_Miles(currFile,path=objPathFolder,returnName = True)
        if type(lamb) == type(None):
            lamb, flux = LD.Load_Dat(currFile,path=objPathFolder) # Dat no da el nombre de la estrella
            if lamb == None:
                print (f"ERROR: FILE {currFile} does not have correct format. Skipping...")
                pass
        if starName != "":
            starName = "_" + starName
        fit,fluxNorm = Normalizar(lamb,flux)    
        LambFluxLists.append((lamb,fluxNorm))
        # Creemos el nombre del fichero
        # Cojamos el nombre del objeto
        name = f"{currFile.split('.')[0]}" + starName
        outNameList.append(name)
    LD.Write_Data(LambFluxLists,outNameList, endPathFolder)
    print(f"\n FILES PROPERLY NORMALISED AND STORED IN \n{endPathFolder}")
    return True