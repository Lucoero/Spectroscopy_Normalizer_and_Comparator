#!/usr/bin/env jupyter
import numpy as np
import scipy as scp
from scipy.signal import medfilt
from scipy.signal import savgol_filter
from scipy.signal import hilbert
from scipy.signal import find_peaks
from scipy.signal import peak_widths
import os as os
import tqdm as tqdm
import matplotlib.pyplot as plt
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
#%% Tipos de Normalizacion
def Norm_Savgol(lamb,flujo, params = [97,4000,-1,"med",1]): #Donde filtro es el tipo de filtro segun filtrar
    """
    Los tipos de filtro posibles son:
        
        med
        mmed
        sg
        h
    params = [proSavgol (97), start, end,filtro,iteraciones]
    El corte es a partir de que longitud de onda (Armstrongs) queremos empezar a normalizar
    ajuste hecho
    flujo normalizado
    """
    try:
        par, start,end,filtro,iteraciones = params
    except:
        print("No se han dado todos los parametros opcionales de Norm_Savgol. Usando parametros por defecto")
        par,start,end,filtro,iteraciones = 97,4000,-1,"med",1
    if end >= lamb[-1]: 
        end = -1
    elif end != -1: 
        end = np.where(lamb > end)[0][0] 
    ajuste = np.copy(flujo)
    start_index = np.where(lamb > start)[0][0]
    
    normalizada_vieja = flujo[start_index:end]
    lamb_cortado = lamb[start_index:end]
    
    for i in range(iteraciones):
        normalizada_nueva = Filtrar(normalizada_vieja,[par], tipo = filtro)
        normalizada_vieja = normalizada_nueva
        
    # Fin ajuste
    ajuste[start_index:end] = normalizada_vieja
    # Nos quitamos el ajuste de 0 porque si no se hacen nans
    ajuste = np.where(ajuste > 0., ajuste, 1.)
    # Ahora normalizamos
    flujo_normalizado = flujo/ajuste
    
    return (lamb_cortado,normalizada_vieja),flujo_normalizado 



def Continuo(fr,pr=0.1,d=5,sg=False,rl=0.5):
    """ # Ya se hace en NormAgg
    if sg:
        f = savgol_filter(fr,5,2)
    else:
        f = fr
    """
    f = fr
    pk1, _ = find_peaks(-1*f,prominence=pr,distance=d)
    w1 = peak_widths(-1*f,pk1,rel_height=rl)
    pk2, _ = find_peaks(f,prominence=pr,distance=d)
    w2 = peak_widths(f,pk2,rel_height=rl)
    pks = np.append(pk1,pk2)
    ws = np.append(w1,w2)
    aux = np.zeros_like(fr)
    for i in range(len(pks)):
        a = int(pks[i] - ws[i])
        b = int(pks[i] + ws[i]) # Porque el pico puede estar en el borde
        
        if b >= len(f)-1: b = len(f)-1
        if b <= 0: b = 0
        if a <= 0: a = 0
        if a >= len(f)-1: a = len(f)-1
        
        lsup = f[b]
        linf = f[a]
        pe = (lsup - linf)/(2*ws[i])
        for t in range(int(a),int(b+1)):
            aux[t] = linf + pe*(t-a)
    for j in range(len(f)):
        if aux[j] == 0:
            aux[j] = f[j]
    return aux,[[pk1,w1],[pk2,w2]]

def Norm_Agg(lamb,flujo,params = [4000,-1,0.1,10,False,0.5]):
    """
    Norm_Agg:
        Normalizamos generando el continuo interpolando por encima de los picos
        los parametros extra son
        start: Donde empezamos el ajuste en Armstrongs
        end: donde terminamos el ajuste en Armstrongs
        pr: Prominencia de los picos
        d: Distancia entre picos
        sg: si usas savgol o no 
        rl: altura relativa
    """
    try:
        start,end,pr,d,sg,rl = params
    except:
        print("No has metido todos los parametros necesarios de Norm_Agg. Usamos los por defecto...")
        start,end,pr,d,sg,rl = 4000,-1,0.1,10,False,0.5
    if sg:
        flujo = savgol_filter(flujo,5,2)
    flujo = np.where(flujo > 0, flujo,0.00001)
    ajuste = np.copy(flujo)
    
    # Buscamos los indices de corte
    start_index = np.where(lamb >= start)[0][0]
    if end >= lamb[-1]: 
        end = -1
    elif end != -1: 
        end = np.where(lamb >= end)[0][0] 
        
    flujo_cortado = flujo[start_index:end]
    lamb_cortado = lamb[start_index:end]
    ajusteLocal = Continuo(flujo_cortado,pr,d,sg,rl)[0] # Sg es usar savgol
    
    # Fin ajuste

    ajuste[start_index:end] = ajusteLocal
    flujo_normalizado = flujo/ajuste
    return (lamb_cortado,ajusteLocal),flujo_normalizado

def Normalise_Folder(objPathFolder,endPathFolder, normFunc, funcParams):
    """
    Normalise_Folder:
        Dado un path con espectros en formato fits o dat,
        y un path de salida, cargamos los datos, los normalizamos,
        y los guardamos en ficheros dat con el nombre original
        y el de la estrella asociado
        Notese que respecto a los fits podemos perder otros datos como la 
        luminosidad, o el corrimiento al rojo, etc.
        
        - normFunc es la funcion de normalizacion que quieres usar (como objeto).
            NORMFUNC DEBE TENER COMO PRIMER PARAMETRO LAMBDA Y COMO SEGUNDO EL FLUJO 
            NORMFUNC DEBE DEVOLVER EN FORMADO 
        - funcParams son los parametros opcionales de la funcion. Si los quieres como defau
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
        ext = os.path.basename(currFile).split(".")[-1]
        if ext == "dat" or ext == "asc":
            lamb, flux = LD.Load_Dat(currFile,path=objPathFolder) # Dat no da el nombre de la estrella
            starName = ""
        elif ext == "fits":
            lamb, flux,starName = LD.Load_Miles(currFile,path=objPathFolder,returnName = True)
        else:
            print (f"ERROR: FILE {currFile} does not have correct format. Skipping...")
            pass
        if starName != "":
            starName = "_" + starName
        fit,fluxNorm = normFunc(lamb,flux, params = funcParams)    
        LambFluxLists.append((lamb,fluxNorm))
        # Creemos el nombre del fichero
        # Cojamos el nombre del objeto
        name = f"{currFile.split('.')[0]}" + starName
        outNameList.append(name)
    LD.Write_Data(LambFluxLists,outNameList, endPathFolder)
    print(f"\n FILES PROPERLY NORMALISED AND STORED IN \n{endPathFolder}")
    return True