# -*- coding: utf-8 -*-
"""
Practica 3
"""
#%% Librerias
import numpy as np

import scipy.optimize as Optimize
import scipy.stats as stats
import scipy.special as special
from scipy.signal import find_peaks
from scipy.signal import peak_widths

import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
#%% Funciones

def Get_Lines_Features(Lamb,Flux,lines,searchWindow = 500):
    """
    Get_Lines_Features:
        Dado un espectro y el diccionario de las lineas deseadas,
        devuelve las propiedades fundamentales de las lineas en forma de 
        [ubicacion,anchura-mitadAltura,profundidad]
    
    TODO: Tener en cuenta tambien lineas de absorcion
    """
    # Recortamos en torno a las lineas para la busqueda
    linesNames = list(lines.keys())
    nLines = len(linesNames)
    linesFeatures = np.zeros((nLines,3)) # [peakPos,peakWidth,peakHeight]
    midWind = searchWindow/2
    lambsCrop = []
    fluxsCrop = []
    for i in range(nLines):
        leftLamb = lines[linesNames[i]]-midWind # Armstrongs
        leftLamb = np.where(Lamb >= leftLamb)[0][0] # Ahora en indices
        rightLamb = lines[linesNames[i]]+midWind
        rightLamb = np.where(Lamb >= rightLamb)[0][0]
        
        # Buscamos los picos en cada corte
        fCrop = Flux[leftLamb:rightLamb]
        lCrop = Lamb[leftLamb:rightLamb]
        peakIndex = find_peaks(-fCrop, prominence = 0.1)[0]
        peakAux = lCrop[peakIndex] # El resto de parametros opcionales, porque deberia ser sencillo encontrarlo

        # Buscamos el pico mas cercano al centro
        indexMin = 0
        minDist = searchWindow
        for j in range(len(peakAux)):
            dist = abs(peakAux[j]-lines[linesNames[i]])
            if dist < minDist:
                minDist = dist
                indexMin = j
        linesFeatures[i,0] = peakAux[indexMin]
        #print(linesFeatures[i,0])
        # Obtenemos anchura
        aux = peak_widths(-fCrop,[peakIndex[indexMin]],rel_height = 0.95)
        width =  aux[0]
        height = aux[1]
        linesFeatures[i,1] = width[0]
        # Obtenemos altura
        linesFeatures[i,2] = height[0]
        # Obtenemos el corte donde se ven estas propiedades
        leftLamb = np.where(Lamb >= linesFeatures[i,0]-linesFeatures[i,1]/2)[0][0]
        rightLamb = np.where(Lamb >= linesFeatures[i,0]+linesFeatures[i,1]/2)[0][0]
        lambsCrop.append(np.array(Lamb[leftLamb:rightLamb]))
        fluxsCrop.append(np.array(Flux[leftLamb:rightLamb]))
    return linesFeatures, lambsCrop,fluxsCrop

def Gauss_Line(x,A,mu,sigma):
    return A*np.exp((x-mu)**2/(2*(sigma)**2))

def Lorentz_Line(x,A,mu,sigma):
    return A*1/(1+(x-mu)**2/(sigma**2))

def Voigt_Line(x,A,mu,gamma): 
    """
    Detalle interesante: Si obtenemos sigma= 0 es lorentziana. gamma = 0, es gaussiana.
    gamma es la anchura a media altura
    """
    return A*special.voigt_profile(x,mu,gamma) 

def Line_Fit(Lamb,Flux, lineCenter,lineWidth,lineHeight):
    """
    Line_Fit:
        Dado un conjunto de datos de flujos y longitudes de onda (que conforman la linea),
        tratamos de buscar una curva de ajuste utilizando una gaussiana con scipy.optimize.curvefit:
            https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.curve_fit.html
        Necesitamos:
            1. El centro de la linea (lo tenemos, idealmente, del diccionario + findpeaks)
            2. Una anchura media (con findpeaks otra vez)
            3. Profundidad de la linea
            4. Una funcion modelo a la que queremos aproximar (en nuestro caso una Gaussiana con 1. y 2,, 
                                                               pero la podemos modificar ligeramente)
    """
    profiles = 3 # Numero de perfiles que usamos
    dArr = np.zeros(profiles) # Array de distancias que obtenemos
    namesArr = [] # Array del nombre de cada metodo
    
    # Aproximamos con scipy a Gaussiana
    try: 
        (gaussA,gaussMu,gaussSigma),gaussErr = Optimize.curve_fit(Gauss_Line,Lamb,Flux,p0 = [lineHeight,lineCenter,lineWidth])
        dArr[0] = stats.wasserstein_distance(Flux, Gauss_Line(Lamb,gaussA,gaussMu,gaussSigma))
        print(gaussA,gaussMu,gaussSigma)
    except RuntimeError: # No converge --> Muy mala aproximacion
        (gaussA,gaussMu,gaussSigma),gaussErr = (0.001,0.001,0.001),np.infty
        dArr[0] = np.infty
    namesArr.append("Gaussian Profile")
    
    # Aproximamos con scipy a Lorentziana
    try:
        (lorA,lorMu,lorSigma),lorErr = Optimize.curve_fit(Lorentz_Line,Lamb,Flux,p0=[lineHeight,lineCenter,lineWidth])
        dArr[1] = stats.wasserstein_distance(Flux, Lorentz_Line(Lamb,lorA,lorMu,lorSigma))
    except RuntimeError:
        (lorA,lorMu,lorSigma),lorErr = (0.001,0.001,0.001),np.infty
        dArr[1] = np.infty
    namesArr.append("Lorentzian Profile")
    
    # Aproximamos con scipy a perfil de Voigt
    try:
        (voigtA,voigtMu,voigtGamma),voigtErr = Optimize.curve_fit(Voigt_Line, Lamb, Flux, p0=[lineHeight,lineCenter,lineWidth]) 
        dArr[2] = stats.wasserstein_distance(Flux, Voigt_Line(Lamb,voigtA,voigtMu,voigtGamma))
    except RuntimeError:
        (voigtA,voigtMu,voigtGamma),voigtErr = (0.001,0.001,0.001),np.infty
        dArr[2] = np.infty
    namesArr.append("Voigt Profile")
    #print(dArr)
    # Tomamos el minimo
    if dArr.all == np.infty:
        print("ERROR: NO METHOD HAS CONVERGED. Returning seedValues ")
        minD = np.infty
        minIndex = -1
        minName = "Unknown"
    else:
        minD = np.min(dArr)
        minIndex = np.argmin(dArr)
        minName = namesArr[minIndex]
        
        print(f"Line_Fit concludes that the best profile is {minName}, with distance D = {minD}")
        print("Returning all profiles data and best index")
    return [[gaussA,gaussMu,gaussSigma], [lorA,lorMu,lorSigma],[voigtA,voigtMu,voigtGamma]], minIndex



#%% Variables entrada 

medB = "B2IV.dat"
bigB = "B2Ib.dat"

lines = {"He I 4922": 4922,r'$H_{\gamma}$': 4340,"Mg II 4481": 4481}

#%% Proceso principal
linesNames = list(lines.keys())
nlines = len(linesNames)
# Cogemos data (ya normalizada)
medLamb, medFlux = LD.Load_Dat(medB,path = "Practica_3")
bigLamb,bigFlux = LD.Load_Dat(bigB, path = "Practica_3")

# Los ploteamos para verlos
SSp.Compare_Spectra([medLamb,bigLamb],[medFlux,bigFlux], NameArr = ["Estrella B9IV (HD027295)", "Estrella B2 Ib (HD206165)"], title = "Estrellas Escogidas Normalizadas", lines = lines)
"""
SSp.Lined_Spectra(bigFlux, bigLamb, lines = lines, title = "Espectro de la estrella HD206165 (B2Ib)")
SSp.Lined_Spectra(medFlux,medLamb,lines = lines, title = "Espectro de la estrella  HD027295 (B9IV)")
"""

# Obtenemos las propiedades de la linea
medFea,medLambsCrop,medFluxsCrop = Get_Lines_Features(medLamb,medFlux,lines)
bigFea,bigLambsCrop,bigFluxsCrop = Get_Lines_Features(bigLamb,bigFlux,lines)

#SSp.Compare_Spectra(medLambsCrop,medFluxsCrop, NameArr = [f"{medFea[0,0]}, teorico: 4922" , f"{medFea[1,0]}, teorico: 4340",f"{medFea[2,0]}, teorico: 4481"])
fits = ["Gaussiana","Lorentziana","Voigt", "No concluyente (no convergencia)"]
# Buscamos el mejor perfil
for i in range(nlines):
    # Obtenemos los parametros para cada 
    print(40*"-" + f"Resultados para la linea {linesNames[i]}:")

    medParams, medIndex = Line_Fit(medLambsCrop[i],medFluxsCrop[i],medFea[i,0],medFea[i,1],medFea[i,2])
    bigParams,bigIndex = Line_Fit(bigLambsCrop[i],bigFluxsCrop[i],bigFea[i,0],bigFea[i,1],bigFea[i,2])
    # Visualizamos
    gaussFlux = Gauss_Line(medLambsCrop[i],medParams[0][0],medParams[0][1],medParams[0][2])
    lorFlux = Lorentz_Line(medLambsCrop[i],medParams[1][0],medParams[1][1],medParams[1][2])
    voigtFlux = Voigt_Line(medLambsCrop[i],medParams[2][0],medParams[2][1],medParams[2][2])
    SSp.Blank_Spectra([medLambsCrop[i],medLambsCrop[i],medLambsCrop[i],medLambsCrop[i]], [medFluxsCrop[i],gaussFlux,lorFlux,voigtFlux],title = f"Comparacion de ajustes de línea {linesNames[i]} para {medB}",multiSpectra = True,spectrasNames = ["Original", "Gaussian","Lorentzian","Voigt"])
    
    gaussFlux = Gauss_Line(bigLambsCrop[i],bigParams[0][0],bigParams[0][1],bigParams[0][2])
    lorFlux = Lorentz_Line(bigLambsCrop[i],bigParams[1][0],bigParams[1][1],bigParams[1][2])
    voigtFlux = Voigt_Line(bigLambsCrop[i],bigParams[2][0],bigParams[2][1],bigParams[2][2])
    SSp.Blank_Spectra([bigLambsCrop[i],bigLambsCrop[i],bigLambsCrop[i],bigLambsCrop[i]], [bigFluxsCrop[i],gaussFlux,lorFlux,voigtFlux],title = f"Comparacion de ajustes de línea {linesNames[i]} para {bigB}",multiSpectra = True,spectrasNames = ["Original", "Gaussian","Lorentzian","Voigt"])
    
    print(f"--> ESTRELLA {medB}, linea {linesNames[i]}:")
    print(f"\t Mejor Fit: {fits[medIndex]}")
    print(f"\t Amplitud (uds): {medParams[i][0]}")
    print(f"\t Posicion (A): {medParams[i][1]}")
    print(f"\t Anchura a media altura (A): {medParams[i][2]}")
    print("")
    print(80*("*"))
    print(f"--> ESTRELLA {bigB}, linea {linesNames[i]}:")
    print(f"\t Mejor Fit: {fits[bigIndex]}")
    print(f"\t Amplitud (uds): {bigParams[i][0]}")
    print(f"\t Posicion (A): {bigParams[i][1]}")
    print(f"\t Anchura a media altura (A): {bigParams[i][2]}")
    print("")