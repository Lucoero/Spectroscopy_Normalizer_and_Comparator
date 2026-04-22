# -*- coding: utf-8 -*-
"""
Script para mostrarnos los espectros de las estrellas de la practica 2.
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as interpolate


import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
import Herramientas as Herr
import normalizar as Norm




plt.close("all")

#%% Tomamos las estrellas

SO= 'O8Ib.fits'
SB= 'B2IV.fits'
SA= 'A5_Ia.fits'
SF= 'F6V.fits'
SG= 'G1_V.fits'

#%% Cargamos los datos como lamb y flux

Lamb1, Flux1 = LD.Load_Miles(SO, path= 'p2EspectrosNormalizar')

Lamb2,Flux2 = LD.Load_Miles(SB, path= 'p2EspectrosNormalizar')

Lamb3,Flux3 = LD.Load_Miles(SA, path= 'p2EspectrosNormalizar')

Lamb4,Flux4 = LD.Load_Miles(SF,path= 'p2EspectrosNormalizar')

Lamb5,Flux5 = LD.Load_Miles(SG,path= 'p2EspectrosNormalizar')

#%% Mostramos las estrellas sin normalizar

Estrella_O= SSp.Blank_Spectra(Lamb1, Flux1, title= 'HD209975, estrella tipo O8 Ia')

Estrella_B= SSp.Blank_Spectra(Lamb2, Flux2, title= 'HD000886, estrella tipo B2 IV')

Estrella_A= SSp.Blank_Spectra(Lamb3, Flux3, title= 'HD017378, estrella tipo A5 Ia')

Estrella_F= SSp.Blank_Spectra(Lamb4, Flux4, title= 'HD206826, estrella tipo F6 V')

Estrella_G= SSp.Blank_Spectra(Lamb5, Flux5, title= 'HD024341, estrella tipo G1 V')





