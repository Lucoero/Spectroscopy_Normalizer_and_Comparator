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




plt.close("all")
#%% Variables de llamada
# Estrellas problema
S1 = "estrella1.dat"
S2 = "estrella2.dat"
S3 = "estrella3.dat"
S4 = "estrella4.dat"

# Espectros del catalogo MILES. https://research.iac.es/proyecto/miles/pages/stellar-libraries/the-catalogue.php
Miles_Name = "s0001.fits"

nCan = 3 # El numero de candidatos que queremos devolver 
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
    r'$H_{\delta}$ (h)': 4100 # Es de las lineas de Fraunhofer. https://arxiv.org/pdf/2410.07301 pg 4
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

MLamb,MFlux = LD.Load_Miles(Miles_Name)
"""
T1 = Par.get_Temp(Lamb1, Flux1)
T2 = Par.get_Temp(Lamb2, Flux2)
T3 = Par.get_Temp(Lamb3, Flux3)
T4 = Par.get_Temp(Lamb4, Flux4)
TArr = np.array([T1,T2,T3,T4])
"""
LambsArr = np.array([Lamb1,Lamb2,Lamb3,Lamb4],dtype = object)
FluxsArr = np.array([Flux1,Flux2,Flux3,Flux4],dtype = object)
#%% Normalizacion
nIter = 50
params = [97] # Es el filtro de mediana (que separacion consideras para cortar y tal)

fit1,N1 = Norm.Normalizar(Lamb1,Flux1,params,iteraciones = nIter)
"""
fit2,N2 = Norm.Normalizar(Lamb2,Flux2,params,iteraciones = nIter)
fit3,N3 = Norm.Normalizar(Lamb3,Flux3,params,iteraciones = nIter)
fit4,N4 = Norm.Normalizar(Lamb4,Flux4,params,iteraciones = nIter)
"""
"""
fitMiles,NMiles = Norm.Normalizar(MLamb,MFlux,iteraciones = 50)
SSp.Compare_Norms(np.array([np.array([MLamb,MFlux],dtype = object)]),np.array([np.array([MLamb,NMiles],dtype = object)]), np.array([fitMiles]))
"""
#%% Busqueda del espectro (KS)
smChosen1,minD1,smCh1,DArr1 = Par.CompareAllSpectra("Catalogo_Miles", (Lamb1,Flux1),lines = lines, distFunc = "WASS", nCandidates = nCan)
#0820, HD208501, B7Iab C
#0711, HD176437, B9.5II-III C
#0718, HD181470, A0III
smChosen2,minD2,smCh2,DArr2 = Par.CompareAllSpectra("Catalogo_Miles", (Lamb2,Flux2),lines = lines,distFunc = "WASS", nCandidates = nCan)
#0817, HD207330, B2.5III
#0252, HD057061, O9II C
#0873, HD219978, K4.5Ib C
smChosen3,minD3,smCh3,DArr3 = Par.CompareAllSpectra("Catalogo_Miles", (Lamb3,Flux3),lines = lines,distFunc = "WASS", nCandidates = nCan)
#0723, CD-24-15398,
#0921, HD107513, Am C
#0473, HD117200, F0 D
smChosen4,minD4,smCh4,DArr4 = Par.CompareAllSpectra("Catalogo_Miles", (Lamb4,Flux4), lines = lines,distFunc = "WASS", nCandidates = nCan)
#0387, HD090508, F9-V C
#0290, HD066573, G5VFe-1.3CH-1 C
#0105, HD018907, G9:V C
#
#%% Ploteado

#SSp.Compare_Spectra(LambsArr,FluxsArr,TArr = TArr,lines = lines)

#SSp.Blank_Spectra(Lamb1,N1)
"""
# Array de todas las normalizaciones para comparar
normal1 = np.array([Lamb1,Flux1],dtype = object)
normal2 = np.array([Lamb2,Flux2],dtype = object)
normal3 = np.array([Lamb3,Flux3],dtype = object)
normal4 = np.array([Lamb4,Flux4],dtype = object)

normalizado1 = np.array([Lamb1,N1],dtype = object)
normalizado2 = np.array([Lamb2,N2],dtype = object)
normalizado3 = np.array([Lamb3,N3],dtype = object)
normalizado4 = np.array([Lamb4,N4],dtype = object)

defArr = np.array([normal1,normal2,normal3,normal4])
normArr = np.array([normalizado1,normalizado2,normalizado3,normalizado4])

fitArr = np.array([fit1,fit2,fit3,fit4],dtype = object)
SSp.Compare_Norms(defArr,normArr,fitArr = fitArr, lines = lines)
"""
'''
# Comparamos los candidatos del primero
normalCompare1 = []
normCompare1 = []
smFit1 = []
normalCompare1.append(np.array([Lamb1,Flux1],dtype = object))
normCompare1.append(np.array([Lamb1,N1],dtype = object))
namesArr1 = []
namesArr1.append("Estrella Problema")
for i in range(len(smCh1)):
    smLamb,smFlux = LD.Load_Miles(smCh1[i])
    normalCompare1.append(np.array([smLamb,smFlux],dtype = object))
    fit, smNorm = Norm.Normalizar(smLamb, smFlux) 
    smFit1.append(fit)
    normCompare1.append(np.array([smLamb,smNorm],dtype = object))
    namesArr1.append(smCh1[i])

normalCompare1 = np.array(normalCompare1)
smFit1 = np.array(smFit1)
normCompare1 = np.array(normCompare1)

SSp.Compare_Norms(normalCompare1,normCompare1, fitArr = smFit1, lines = lines, title = "Comparacion para el espectro problema 1",TArr= namesArr1)
'''
