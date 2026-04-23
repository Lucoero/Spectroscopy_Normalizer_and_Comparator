# -*- coding: utf-8 -*-
"""
Main
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



plt.close("all")
#%% Variables de llamada
# Estrellas problema
S1 = "estrella1.dat"
S2 = "estrella2.dat"
S3 = "estrella3.dat"
S4 = "estrella4.dat"
problemNStars = "Estrellas_Problema" # Carpeta donde estan las normalizaciones
SN1Name = "SN1.dat"
SN2Name = "SN2.dat"
SN3Name = "SN3.dat"
SN4Name = "SN4.dat"

nCan = 3 # El numero de candidatos que queremos devolver con la funcion comparador

useAgg = True # Si queremos usar la normalizacion agresiva o no 
# El diccionario de lineas espectrales esta en LinesLib.py
lines = LinesLib.lines
#%% Obtencion de datos

Lamb1, Flux1 = LD.Load_Dat(S1)

Lamb2,Flux2 = LD.Load_Dat(S2)

Lamb3,Flux3 = LD.Load_Dat(S3)

Lamb4,Flux4 = LD.Load_Dat(S4)


#MLamb,MFlux = LD.Load_Miles(Miles_Name)
LambsList = [Lamb1,Lamb2,Lamb3,Lamb4]
FluxsList = [Flux1,Flux2,Flux3,Flux4]
#%% Normalizacion
# Llamamos a las normalizaciones que tenemos
Lamb1,SN1 = LD.Load_Dat(SN1Name,problemNStars)
Lamb2,SN2 = LD.Load_Dat(SN2Name,problemNStars)
Lamb3,SN3 = LD.Load_Dat(SN3Name,problemNStars)
Lamb4,SN4 = LD.Load_Dat(SN4Name,problemNStars)

NormProblema = [(Lamb1,SN1),(Lamb2,SN2),(Lamb3,SN3),(Lamb4,SN4)]

#%% Busqueda del espectro

smChosen1,minD1,smCh1,DArr1 = Par.CompareAllSpectra("IgnacioNormalizado", (Lamb1,SN1),lines = lines, distFunc = "WASS", nCandidates = nCan)
# CON VIEJA NORMALIZACION
#0820, HD208501, B7Iab C
#0711, HD176437, B9.5II-III C
#0718, HD181470, A0III

# CON NUEVA NORMALIZACION + Catalogo Ignacio
# HD24398 --> B2 o B4
# El problema es que lo ha cogido porque tiene un intervalo para la norma mucho menor

smChosen2,minD2,smCh2,DArr2 = Par.CompareAllSpectra("IgnacioNormalizado", (Lamb2,SN2),lines = lines,distFunc = "WASS", nCandidates = nCan)
# CON VIEJA NORMALIZACION
#0817, HD207330, B2.5III
#0252, HD057061, O9II C
#0873, HD219978, K4.5Ib C

# CON NUEVA NORMALIZACION + Catalogo Ignacio
#HD175191, B2.5V
#HD148605, B3V
# HD36936 B3/5V
# HD033276 A1V
# HD027524 ??
smChosen3,minD3,smCh3,DArr3 = Par.CompareAllSpectra("MilesNormalizado", (Lamb3,SN3),lines = lines,distFunc = "WASS", nCandidates = nCan)
# CON VIEJA NORMALIZACION
#0723, CD-24-15398, K0??
#0921, HD107513, Am C
#0473, HD117200, F0 D
# CON NUEVA NORMALIZACION
#HD137391; F0IV / G1V

smChosen4,minD4,smCh4,DArr4 = Par.CompareAllSpectra("MilesNormalizado", (Lamb4,SN4), lines = lines,distFunc = "WASS", nCandidates = nCan)
# CON VIEJA NORMALIZACION
#0387, HD090508, F9-V C
#0290, HD066573, G5VFe-1.3CH-1 C
#0105, HD018907, G9:V C
# CON NUEVA NORMALIZACION
# HD024451, K4
# HD036003, K5V
# HD029065, K4 III
#%% Ploteado

#SSp.Compare_Spectra(LambsList,FluxsList,NameArr = [S1,S2,S3,S4],lines = lines, title = "Estrellas Problema")

#%%% Ploteado de los resultados de la busqueda
print("Veamos ahora todos los espectros de todos los candidatos")

# Comparamos los candidatos del primero
lambCompare1 = []
fluxCompare1 = []

lambCompare1.append(Lamb1)
fluxCompare1.append(SN1)
namesArr1 = []
namesArr1.append(S1)
for i in range(len(smCh1)):
    smLamb,smFlux = LD.Load_Dat(smCh1[i],path = "IgnacioNormalizado")
    lambCompare1.append(smLamb) 
    fluxCompare1.append(smFlux)
    namesArr1.append(smCh1[i])
SSp.Compare_Spectra(lambCompare1 ,fluxCompare1, lines = lines, title = "Comparacion para el espectro problema 1",NameArr= namesArr1)


# Comparamos los candidatos del segundo

lambCompare2 = []
fluxCompare2 = []

lambCompare2.append(Lamb2)
fluxCompare2.append(SN2)

namesArr2 = []
namesArr2.append(S2)
for i in range(len(smCh2)):
    smLamb,smFlux = LD.Load_Dat(smCh2[i], path = "IgnacioNormalizado")
    lambCompare2.append(smLamb) 
    fluxCompare2.append(smFlux)
    namesArr2.append(smCh2[i])
SSp.Compare_Spectra(lambCompare2 ,fluxCompare2, lines = lines, title = "Comparacion para el espectro problema 2",NameArr= namesArr2)

# Comparamos los candidatos del Tercero

lambCompare3 = []
fluxCompare3 = []

lambCompare3.append(Lamb3)
fluxCompare3.append(SN3)

namesArr3 = []
namesArr3.append(S3)
for i in range(len(smCh3)):
    smLamb,smFlux = LD.Load_Dat(smCh3[i], path = "MilesNormalizado")
    lambCompare3.append(smLamb) 
    fluxCompare3.append(smFlux)
    namesArr3.append(smCh3[i])
SSp.Compare_Spectra(lambCompare3 ,fluxCompare3, lines = lines, title = "Comparacion para el espectro problema 3",NameArr= namesArr3)

# Comparamos los candidatos del Cuarto

lambCompare4 = []
fluxCompare4 = []

lambCompare4.append(Lamb4)
fluxCompare4.append(Flux4)

namesArr4 = []
namesArr4.append(SN4)
for i in range(len(smCh4)):
    smLamb,smFlux = LD.Load_Dat(smCh4[i], path = "MilesNormalizado")
    lambCompare4.append(smLamb) 
    fluxCompare4.append(smFlux)
    namesArr4.append(smCh4[i])
SSp.Compare_Spectra(lambCompare4 ,fluxCompare4, lines = lines, title = "Comparacion para el espectro problema 4",NameArr= namesArr4)
