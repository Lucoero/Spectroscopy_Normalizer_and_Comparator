# -*- coding: utf-8 -*-
"""
NormalizarEstrella Concreta:
    Coges una estrella en concreto, cambiamos el parametro de altura en lectura de pico,
    y se normaliza y se guarda
"""
import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
import Herramientas as Herr
import normalizar as Norm
import LinesLib as LinesLib
#%% Variables entrada
isDat = False
starFile = "B2Ib.fits"
starDir = "Practica_3"
outDir = "Practica_3"
outFilename = "B2Ib"

figTitle = "Estrella Tipo B2 IV"
rl = 0.85 # Altura relativa de los picos empezando desde el fondo desde el cual mides la altura. LO BASICO A CAMBIAR

startNorm = 3930 # Lamb donde quieres empezar la normalizacion
endNorm = 7000 # Indice donde quieres acabar la normalizacion
nIter = 50 # N iteraciones si usas filtro savgol
savgolParams = [97,startNorm,endNorm, "med",nIter]
"""
Agg Params:
    startNorm: El lambda en el que empezamos el ajuste (A)
    endNorm: El lambda en el que acabamos el ajuste (A). Si se usa -1 llega hasta el final
    pr: Prominencia de lo que se considera pico
    d: distancia entre picos
    sg: Si queremos usar filtro savgol para eliminar ruido antes de normalizar
    rl: Altura relativa de los picos empezando desde el fundo del pico
"""
aggParams = [startNorm,endNorm,0.08,10,False,rl] # El ultimo parametro es la altura relativa donde interpolamos el pico (de abajo a arriba pico. 0 Es coger el fondo del pico)
useAgg = True

# Si quieres normalizar una carpeta entera
normFolder = "Catalogo_Ignacio"
outputFolder = "IgnacioNormalizado"
#%% Proceso principal
lines = LinesLib.lines

# Cargo info
if isDat:
    lamb,flux = LD.Load_Dat(starFile,path = starDir)
else:
    lamb,flux = LD.Load_Miles(starFile,path = starDir)

# Normalizo
if useAgg:
    fit, fluxN = Norm.Norm_Agg(lamb,flux,aggParams)
else:
    fit,fluxN = Norm.Norm_Savgol(lamb,flux,aggParams)

# Veo la normalizacion
normSpectra = [lamb,fluxN]
originalSpectra = [lamb,flux]

SSp.Compare_Norms([originalSpectra],[normSpectra],fitArr = [fit],NameArr = ["S1"], title = figTitle)# lines = lines)
# Guardo
LD.Write_Data([normSpectra],[outFilename],outDir)

#%% Si queremos normalizar un folder entero
#aggParams[-1] = 0.5 # Prominencia siempre 0.5
Norm.Normalise_Folder(normFolder, outputFolder, Norm.Norm_Agg, aggParams)