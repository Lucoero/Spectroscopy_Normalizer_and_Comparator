# -*- coding: utf-8 -*-
"""
Practica 3
"""
#%% Librerias
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as Interpolate
import scipy.optimize as Optimize
import scipy.stats as stats
import scipy.special as special

import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
import Herramientas as Herr
import normalizar as Norm
#%% Variables entrada

medB = "B2IV.dat"
bigB = "B2Ib.dat"

lines = {"He I 4922": 4922,r'$H_{\gamma}$': 4340,"Mg II 4481": 4481}


#%% Proceso principal

medFlux, medLamb = LD.Load_Dat(medB,path = "Practica_3")
bigFlux,bigLamb = LD.Load_Dat(bigB, path = "Practica_3")

# Los ploteamos para verlos
SSp.Compare_Spectra([medFlux,bigFlux],[medLamb,bigLamb], NameArr = ["Estrella B9IV (HD027295)", "Estrella B2 Ib (HD206165)"], title = "Estrellas Escogidas Normalizadas", lines = lines)
"""
SSp.Lined_Spectra(bigFlux, bigLamb, lines = lines, title = "Espectro de la estrella HD206165 (B2Ib)")
SSp.Lined_Spectra(medFlux,medLamb,lines = lines, title = "Espectro de la estrella  HD027295 (B9IV)")
"""
# Normalizamos

# Vemos normalizaciones

# Obtenemos las propiedades de la linea

# Aislamos las lineas que queremos y las aproximamos
lambMet = lines["Ca II (K)"]
lambHe = lines["He I 4922"]
lambBa = lines[r'$H_{\gamma}$']


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
    namesArr = np.zeros(profiles,dtype = str) # Array del nombre de cada metodo
    
    # Aproximamos con scipy a Gaussiana
    (gaussA,gaussMu,gaussSigma),gaussErr = Optimize.curve_fit(Gauss_Line,Lamb,Flux,[lineHeight,lineCenter,lineWidth])
    dArr[0] = stats.wasserstein_distance(Flux, Gauss_Line(Lamb,gaussA,gaussMu,gaussSigma))
    namesArr[0] = "Gaussian Profile"
    
    # Aproximamos con scipy a Lorentziana
    (lorA,lorMu,lorSigma),lorErr = Optimize.curve_fit(Lorentz_Line,Lamb,Flux,[lineHeight,lineCenter,lineWidth])
    dArr[1] = stats.wasserstein_distance(Flux, Lorentz_Line(Lamb,lorA,lorMu,lorSigma))
    namesArr[1] = "Lorentzian Profile"
    
    # Aproximamos con scipy a perfil de Voigt
    voigtA,voigtMu,voigtGamma,voigtErr = Optimize.curve_fit(Voigt_Line, Lamb, Flux, [lineHeight,lineCenter,lineWidth]) 
    dArr[2] = stats.wasserstein_distance(Flux, Voigt_Line(Lamb,voigtA,voigtMu,voigtGamma))
    namesArr[2] = "Voigt Profile"
    
    # Tomamos el minimo
    minD = np.min(dArr)
    minIndex = np.where(dArr == minD)
    minName = namesArr[minIndex]
    
    # Visualizamos
    gaussFlux = Gauss_Line(Lamb,gaussA,gaussMu,gaussSigma)
    lorFlux = Lorentz_Line(Lamb, lorA, lorMu, lorSigma)
    voigtFlux = Voigt_Line(Lamb, voigtA, voigtMu, voigtGamma)
    SSp.Blank_Spectra([Lamb,Lamb,Lamb,Lamb], [Flux,gaussFlux,lorFlux,voigtFlux],title = "Distintos perfiles de líenas obtenidos",multiSpectra = True,spectrasNames = ["Original", "Gaussian","Lorentzian","Voigt"])
    
    print(f"Line_Fit concludes that the best profile is {minName}, with distance D = {minD}")
    print("Returning all profiles data and best index")
    return [[gaussA,gaussMu,gaussSigma], [lorA,lorMu,lorSigma],[voigtA,voigtMu,voigtGamma]], minIndex


