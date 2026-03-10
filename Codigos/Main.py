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
Lineas= {
    'Hb': 6563,
    'Hgamma': 4861,
    'Hdelta': 4340,
    'Hepsilon': 4120
    }

#%% Obtencion de datos

Lamb1, Flux1 = LD.Load_Dat(S1)

Lamb2,Flux2 = LD.Load_Dat(S2)
#MLamb,MFlux = LD.Load_Miles(Miles_Name)
#%% Correcciones
#%% Ploteado

SSp.Compare_Spectra(Lamb1,Flux1,Lamb2,Flux2)

#SSp.Blank_Spectra(MLamb, MFlux)