# -*- coding: utf-8 -*-
"""
Main
"""
import numpy as np
import matplotlib.pyplot as plt
import Load_Data as LD
import Show_Spectra as SSp

plt.close("all")
#%% Variables de llamada
# Estrellas problema
S1 = "estrella1.dat"
S2 = "estrella2.dat"
S3 = "estrella3.dat"
S4 = "estrella4.dat"

# Espectros del catalogo MILES
Miles_Name = "s0010.fits"


# Diccionario de lineas que vamos a marcar
lineas= {
    'Hb': 6563,
    'Hgamma': 4861,
    'Hdelta': 4340,
    'Hepsilon': 4120
    }

#%% Obtencion de datos

Lamb1, Flux1 = LD.Load_Dat(S1)

Lamb2,Flux2 = LD.Load_Dat(S2)

Lamb3,Flux3 = LD.Load_Dat(S3)

Lamb4,Flux4 = LD.Load_Dat(S4)

#MLamb,MFlux = LD.Load_Miles(Miles_Name)

T1 = 

#%% Correcciones
#%% Ploteado

SSp.Compare_Spectra(np.array([Lamb1,Lamb2,Lamb3,Lamb4],dtype = object),np.array([Flux1,Flux2,Flux3,Flux4],dtype = object),lines = lineas)

#SSp.Blank_Spectra(MLamb, MFlux)