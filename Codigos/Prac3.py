# -*- coding: utf-8 -*-
"""
Practica 3
"""
#%% Librerias
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as Interpolate
import scipy.optimize as Optimize

import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
import Herramientas as Herr
import normalizar as Norm
#%% Variables entrada

medB = "B9IV.fits"
bigB = "B2Ib.fits"


lineas_metalicas= {
    "Ca II (K)": 3934,
    }

lineas_balmer= {
    r'$H_{\gamma}$': 4340,
    }

lineas_helio= {
    "He I 4922": 4922,
    }

lines = {"He I 4922": 4922,r'$H_{\gamma}$': 4340,"Ca II (K)": 3934}


#%% Proceso principal

medFlux, medLamb = LD.Load_Miles(medB,path = "Practica_3")
bigFlux,bigLamb = LD.Load_Miles(bigB, path = "Practica_3")

# Los ploteamos para verlos
SSp.Compare_Spectra([medFlux,bigFlux],[medLamb,bigLamb], NameArr = ["Estrella B9IV (HD027295)", "Estrella B2 Ib (HD206165)"], title = "Estrellas Escogidas", lines = lines)
"""
SSp.Lined_Spectra(bigFlux, bigLamb, lines = lines, title = "Espectro de la estrella HD206165 (B2Ib)")
SSp.Lined_Spectra(medFlux,medLamb,lines = lines, title = "Espectro de la estrella  HD027295 (B9IV)")
"""
# Normalizamos

# Vemos normalizaciones

# Aislamos las lineas que queremos y las aproximamos
def Gauss_Line(x,A,mu,sigma):
    return A*np.exp((x-mu)**2/(2*(sigma)**2))
    
def Line_Fit(Lamb,Flux):
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
    # Buscamos el centro de la linea, su anchura media y profundidad como semilla base
    mu = 0 # Centro linea
    sigma = 0 # Anchura linea
    A = 0 # Altura linea
    # Aproximamos con scipy
    bestA,bestMu,bestSigma = Optimize.curve_fit(Gauss_Line,Lamb,Flux,[A,mu,sigma])[0]
    
    # Visualizamos
    bestFlux = Gauss_Line(Lamb,bestA,bestMu,bestSigma)
    SSp.Blank_Spectra([Lamb,Lamb], [Flux,bestFlux],title = "Ajuste obtenido",multiSpectra = True,spectrasNames = ["Original", "Ajuste"])
    return bestA,bestMu,bestSigma


