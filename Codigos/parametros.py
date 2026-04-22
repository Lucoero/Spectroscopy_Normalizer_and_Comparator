"""
parametros:
    Se encarga de tomar los archivos y calcular
    diversos parametros como la Temperatura y demas
    
Alguna libreria de parametros util
https://sites.uni.edu/morgans/astro/course/Notes/section2/spectraltemps.html

"""

import os as os
import re
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.visualization import quantity_support
from pathlib import Path
from scipy.signal import convolve
from scipy.signal import find_peaks, peak_prominences, peak_widths
import tqdm as tqdm

from scipy.interpolate import make_smoothing_spline
from scipy.stats import kstest
import scipy.spatial.distance as distances
import scipy.stats as stats
import Load_Data as Load
import normalizar as normalizar
import Show_Spectra as SSp
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
#%% Calculo de picos
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
#%% Comparacion automatica
def CompareAllSpectra(dataFolder,objSpectra, outFolder = "Outputs", lines = {}, distFunc = "WASS", nCandidates = 1):
    """
    CompareAllSpectra
    Dada una carpeta con ficheros en formato miles ([lambda,flux]), y un espectro problema en formato
    ([lambda, flux]), ambos normalizados. Se realiza el siguiente proceso:
        1. Interpolamos total (con picos y todo) el espectro problema (sp).
        2. Cogemos un espectro de miles (sm). Lo normalizamos y le realizamos una interpolacion total (con picos y todo)
        3. Calculamos la distancia entre las distribuciones. Usamos un dominio con un numero total de puntos nPoints Almacenamos el valor de la distancia en un array [Di]
        4. Cogemos el siguiente espectro y repetimos (2,3)
        6. Cogemos el minimo del array. Buscamos su sm correspondiente por indice
        
    Si KS no va bien podemos usar Chi^2
    Renormalizamos los espectros
    """
    #%%% Obtencion de datos
    print(60*"*" + "\nPREPARING COMPARISON WITH SPECTRAS" +"\n" +60*"*")
    spLamb = objSpectra[0].copy()
    spFlux = objSpectra[1].copy()
    if not os.path.exists(dataFolder):
        print(f"ERROR: No existe el directorio con archivos: \n {dataFolder}")
        return None
    FilesArr = os.listdir(dataFolder)
    n = len(FilesArr)
    DArr = np.zeros(n)
    #pArr = np.zeros(n) # Array de p-values relacionados con el test.
    #%%% Ajuste del espectro problema
    # Ajusto el objetivo con splines
    """
    ASUNTO IMPORTANTE: SI EL DOMINIO DE LAMBDAS DEL MILES ES MENOR AL DE LOS ESPECTROS
    OBJETIVO, AL HACER EL AJUSTE ESTE DEBE EXTRAPOLAR EL MILES HASTA LAS LAMBDAS NECESARIAS.
    ESTO HACE QUE HAYA UNA "COLA" EXPONENCIAL QUE METE ERROR.
    ASI QUE LO QUE HAY QUE HACER ES TRUNCAR LAS LAMBDAS, NUNCA EXTRAPOLAR
    """
    #spInterpol = make_smoothing_spline(spLamb, spFlux,lam = 0)
    spMinLamb = np.min(spLamb)
    spMaxLamb = np.max(spLamb)
    #%%% Comparando con la libreria de espectros
    # Voy archivo a archivo de Miles analizando
    for i in tqdm.tqdm(range(n)):
        currFile = FilesArr[i]
        extension = currFile.split(".")[-1]
        if extension == "dat" or extension == "asc":
            smLamb,smFlux = Load.Load_Dat(currFile, path = dataFolder)
        elif extension == "fits":   
            smLamb,smFlux = Load.Load_Miles(currFile, path = dataFolder)
        else:
            print(f"ERROR: LA EXTENSION NO CORRESPONDE A NINGUNA COMPATIBLE: {extension}")
            return None
        # Compruebo si el espectro problema o el objetivo tiene menor rango de lambdas
        minLamb = min(np.min(smLamb),spMinLamb)
        maxLamb = max(np.max(smLamb),spMaxLamb)
        minIndex = np.where(spLamb > minLamb)[0][0]
        maxIndex = np.where(spLamb >= maxLamb)[0][0]
        lambArr = spLamb[minIndex:maxIndex]
        # Ajusto con splines
        smInterpol =  make_smoothing_spline(smLamb, smFlux,lam = 0)
        # Comparo  en un mismo dominio
        
        if distFunc == "KS":
            result = kstest(spFlux,smInterpol(lambArr))[0] # Podriamos guardar tambien los pvalues
        elif distFunc == "WASS":
            result = stats.wasserstein_distance(spFlux, smInterpol(lambArr))
        else:
            print("No coincide con ninguna")
            result = 0
        # Almaceno los resultados
        DArr[i] = result 
        
    # Organizo el array de minimos para ordenarlos de menor a mayor
    DSorted = np.sort(DArr)
    minD = DSorted[:nCandidates] # Tomamos N candidatos con la minima distancia posible
    smChosen = []
    for i in range(nCandidates):
        currIndex = np.where(DArr == DSorted[i])[0][0]
        smChosen.append(FilesArr[currIndex])
    smChosen = np.array(smChosen)
    # Printeo el resultado
    print(f"\nTests says that {smChosen[0]} is the most probable spectra with distance {minD[0]}. {distFunc} was used \n The {nCandidates} Best Candidates are in the following order:")
    for i in range(nCandidates):
        print(f"{i+1}: {smChosen[i]};  D = {minD[i]}") 
    print("Showing the comparison of spectras for the minimun distance")
    # Printeamos la comparacion
    currFile = smChosen[0]
    extension = currFile.split(".")[-1]
    if extension == "dat" or extension  == "asc":
        smLamb,smFlux = Load.Load_Dat(currFile, path = dataFolder)
    elif extension == "fits": 
        smLamb,smFlux = Load.Load_Miles(currFile, path = dataFolder)
    SSp.Compare_Spectra([objSpectra[0],smLamb], [objSpectra[1],smFlux],title = f"Candidato para espectro: {smChosen[0]} con D = {minD[0]}", lines=lines)
    return smChosen,minD,smChosen,DArr