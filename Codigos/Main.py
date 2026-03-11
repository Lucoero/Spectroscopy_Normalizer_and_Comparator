# -*- coding: utf-8 -*-
"""
Main
"""
import numpy as np
import matplotlib.pyplot as plt

import Load_Data as LD
import Show_Spectra as SSp
import parametros as Par
import Herramientas as Herr
import normalizar as norm




plt.close("all")
#%% Variables de llamada
# Estrellas problema
S1 = "estrella1.dat"
S2 = "estrella2.dat"
S3 = "estrella3.dat"
S4 = "estrella4.dat"

# Espectros del catalogo MILES. https://research.iac.es/proyecto/miles/pages/stellar-libraries/the-catalogue.php
Miles_Name = "s0010.fits"


# Diccionario de lineas que vamos a marcar
"""
Nota: Lineas mas destacables

    Balmer(partes segundo nivel) https://es.wikipedia.org/wiki/Líneas_de_Balmer
        - Halpha (salto 2 --> 3): 6563 A
        -Hbeta (salto 2 --> 4): 4861 
        -HY (2--> 5): 4340 A
        -Hdelta (2-->6): 4107
        -Hepsi (2--> 7): 3970
        -Hchi (2-->8): 3889
        -Heta (2--> 9): 3835
        
    Lymann(partes del primer nivel) no las vemos en el optico, por eso no las tenemos en cuenta 
    (son de 1000 A)
    
    Calcio H y K 
    
    Helio I y He II
    
    Otros metales (Fe,Mg,Si)
"""

lineas_metalicas= {
    "Ca I 4227": 4227,
    "Ca II (K)": 3934,
    "Ca II (H)": 3968,
    "Mg II 4481": 4481,
    "Si IV 4654": 4654,
    "Si IV 4631": 4631,
    "Ti II 4179":4179,
    "Na (D2)": 5890
    }

lineas_fe= {
    "Fe I 4383": 4383,
    "Fe I 4271": 4271,
    "Fe II 4233": 4233,
    "Fe II 4172": 4172
    }

lineas_balmer= {
    r'$H_{\alpha}$': 6563,
    r'$H_{\gamma}$': 4861,
    r'$H_{\delta}$': 4340,
    r'$H_{\epsilon}$': 4120,   
    }

lineas_helio= {
    "He I 5875": 5875,
    "He I 4922": 4922,
    "He I 4713": 4713,
    "He I 4471": 4471,
    "He I 4388": 4388,
    "He I 4144": 4144,
    "He II 4686": 4686,
    "He II 4542": 4542,
    "He II 4200": 4200,
    "He I + II 4026": 4026  
    }

lines= {}
lines.update(lineas_metalicas)
lines.update(lineas_balmer)
lines.update(lineas_helio)
lines.update(lineas_fe)

#%% Obtencion de datos

Lamb1, Flux1 = LD.Load_Dat(S1)

Lamb2,Flux2 = LD.Load_Dat(S2)

Lamb3,Flux3 = LD.Load_Dat(S3)

Lamb4,Flux4 = LD.Load_Dat(S4)

#MLamb,MFlux = LD.Load_Miles(Miles_Name)

T1 = Par.get_Temp(Lamb1, Flux1)
T2 = Par.get_Temp(Lamb2, Flux2)
T3 = Par.get_Temp(Lamb3, Flux3)
T4 = Par.get_Temp(Lamb4, Flux4)
TArr = np.array([T1,T2,T3,T4])

LambsArr = np.array([Lamb1,Lamb2,Lamb3,Lamb4],dtype = object)
FluxsArr = np.array([Flux1,Flux2,Flux3,Flux4],dtype = object)
#%% Normalizacion
'''
start= np.where(Flux1 == np.max(Flux1))
Nada, Lamb1_lineas,Nada2= Herr.BuscadorMinimos(np.array([Lamb1,Flux1],dtype= object),start)
#%% Ploteado
'''

SSp.Compare_Spectra(LambsArr,FluxsArr,TArr = TArr,lines = lines)
plt.show()

#SSp.Blank_Spectra(MLamb, MFlux)
par = [93]
y = norm.filtro(3,Flux1, par)


plt.plot(Lamb1, Flux1)
plt.plot(Lamb1, y)
plt.show()
plt.plot(Lamb1,[Flux1[i]/y[i] for i in range(len(Flux1))])
