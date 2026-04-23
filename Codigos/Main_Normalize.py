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
starFile = "F6V.fits"
starDir = "p2EspectrosNormalizar"
outDir = "p2EspectrosNormalizar"
outFilename = "Estrella_F"

figTitle = "HD206826, estrella tipo F6 V Normalizada"
rl = 0.85 # Altura relativa de los picos empezando desde el fondo desde el cual mides la altura. LO BASICO A CAMBIAR

startNorm = 4000 # Lamb donde quieres empezar la normalizacion
endNorm = 7500 # Indice donde quieres acabar la normalizacion
nIter = 50 # N iteraciones si usas filtro savgol
savgolParams = [97,startNorm,endNorm, "med",nIter]
"""
Agg Params:
    startNorm: El lambda en el que empezamos el ajuste (A)
    endNorm: El lambda en el que acabamos el ajuste (A). Si se usa -1 llega hasta el final
    pr: Prominencia de lo que se considera pico
    sg: Si queremos usar filtro savgol para eliminar ruido antes de normalizar
    rl: Altura relativa de los picos empezando desde el fundo del pico
"""
aggParams = [startNorm,endNorm,0.1,10,False,rl] # El ultimo parametro es la altura relativa donde interpolamos el pico (de abajo a arriba pico. 0 Es coger el fondo del pico)
useAgg = True

# Si quieres normalizar una carpeta entera
normFolder = "Catalogo_Miles"
outputFolder = "MilesNormalizado"
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

SSp.Compare_Norms([originalSpectra],[normSpectra],fitArr = [fit],NameArr = ["S1"], title = figTitle) #lines = lines)
# Guardo
LD.Write_Data([normSpectra],[outFilename],outDir)

#%% Si queremos normalizar un folder entero
#aggParams[-1] = 0.5 # Prominencia siempre 0.5
#Norm.Normalise_Folder(normFolder, outputFolder, Norm.Norm_Agg, aggParams)