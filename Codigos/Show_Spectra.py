# -*- coding: utf-8 -*-
"""
Show_Spectra

Dado los arrays de datos, ploteamos los espectros con distintas propiedades
"""

import numpy as np
import matplotlib.pyplot as plt

#%% Variables globales de plots
yscale = 1.01 # Escala para escoger el tamanno de los plots
lWidth = 1 # Grosor de las lineas de plot
lineScale = 1 # Escala para las lineas atomicas
lineSize = 8 # Tamaño del texto de las lineas
ncol = 1 # Numero de columnas para la leyenda

fontsize = 11 #  text scale for plots
numbersize = [13,13] # scale for number axes
plt.rc('font', size=fontsize) # Para el escalado de texto en plots
plt.xticks(fontsize = numbersize[0])
plt.yticks(fontsize = numbersize[1])
plt.close("all")
#%% Funciones
def Pad_Array(arr): # Gracias, https://stackoverflow.com/questions/24494356/how-to-find-min-max-values-in-array-of-variable-length-arrays-with-numpy   
    nArr = len(arr)
    M = max(len(a) for a in arr) # Calculamos la longitud maxima de los arrays
    out = np.zeros((nArr,M),dtype = object)
    for i in range(nArr):
        out[i] = np.array(list(arr[i])+[np.nan]*(M-len(arr[i])))
    return out


def Blank_Spectra(lambList,fluxList, title = "Espectro", multiSpectra = False,spectrasNames = []):
    """
    Blank_Spectra:
        
        Ploteamos el espectro con la menor cantidad de cambios posible.
        Perfecto para un primer acercamiento
    """
    fig, ax = plt.subplots(figsize = (15,6))
    fig.suptitle(title)
    Axe_Blank_Spectra(lambList, fluxList, ax, multiSpectra = multiSpectra, spectrasNames = spectrasNames)
    fig.show()
    return

def Axe_Blank_Spectra(lambList,fluxList,ax, name = False, show_yName = True, show_xName = True, multiSpectra = False, spectrasNames = []):
    if show_xName:
        ax.set_xlabel(r"$\lambda \ (\mathring{A})$")
    if show_yName:
        ax.set_ylabel("Flux (uds)")
    if multiSpectra: # Tenemos un array de lambdas de distintos espectros que queremos superponer
        nSpectra = len(lambList)
        if spectrasNames == []:
            spectrasNames = [f"{i}" for i in range(nSpectra)]
        for i in nSpectra:
            ax.plot(lambList[i],fluxList[i], label = spectrasNames[i])
            ax.legend()
    ax.plot(lambList,fluxList, linewidth = lWidth)
    if name:
        ax.set_title(name,loc="right", y=.5,
        rotation=270, ha="left", va="center")
    return 
def Lined_Spectra(lamb,flux,lines, title = "Lined_Spectra", show_yName = True,multiSpectra = False, spectrasNames = []):
    """
    Lined_Spectra:
        Ploteamos el espectro con las lineas que queramos marcar.
        Lines debe enviarse como un diccionario con {"line_name" : lamb,}
    """
    # Ploteamos espectro
    fig,ax = plt.subplots(figsize = (15,6))
    fig.suptitle(title)
    Axe_Lined_Spectra(lamb,flux,lines,ax, show_yName = show_yName,multiSpectra = multiSpectra, spectrasNames = spectrasNames)
    fig.legend(ncol = ncol)
    fig.show()
    return
 
def Axe_Lined_Spectra(lamb,flux,lines,ax, name = False, show_yName = True, show_xName = True, show_LineName = True, multiSpectra = False, spectrasNames = []):
    Axe_Blank_Spectra(lamb, flux, ax, name = name, show_yName = show_yName, show_xName = show_xName, multiSpectra = multiSpectra, spectrasNames=spectrasNames)
    minLine = np.min(flux)*yscale
    maxLine = np.max(flux)*yscale
    # Ploteamos las lineas
    for name in lines:
        if show_LineName:
            ax.plot((lines[name],lines[name]), (minLine,maxLine),label = name,linestyle = "dashed",linewidth = lWidth * lineScale) # Este para el label
            ax.annotate(name,xy = (lines[name] + 5,maxLine -0.1*maxLine), xycoords = "data", ha = 'center', va = 'bottom',rotation = 'vertical', size = lineSize)
        else:
            ax.plot((lines[name],lines[name]), (minLine,maxLine),linestyle = "dashed",linewidth = lWidth * lineScale) # Este para el label
    return 
    
def Compare_Spectra(lambArr,fluxArr,NameArr = [],lines = {}, title = "Comparison between Spectras", show_yName = True, saveFig = False):
    """
    Compare_Spectra:
        Ploteamos uno al lado del otro varios espectros para compararlos. Si se dan lineas
        como parametro opcional, las coloca tambien en ambos
        
    Entrada:
        LambsArr: Array de arrays de longitudes de onda
        FluxArr: Array de arrays de Flujos
        TArr: Array de temperaturas de cada una
    """
    n = len(lambArr)
    
    fig,ax = plt.subplots(n,1, figsize = (60,6),sharex=True, sharey=False)
    fig.subplots_adjust(hspace=0.001)
    fig.suptitle(title)
    Axe_Compare_Spectra(lambArr,fluxArr,ax, NameArr = NameArr, lines = lines, show_yName = show_yName)
    fig.legend(ncol = ncol)
    if saveFig:
        fig.savefig("Figuras/" + title + ".jpg",format='jpg')
    fig.show()
    return

def Axe_Compare_Spectra(lambArr,fluxArr, ax,NameArr = [],lines = {}, show_yName = False):
    """
    Axe_Compare_Spectra:
        Hace lo que hace Compare_Spectra, pero para utilizarse por otros plots
    """
    n = len(ax)
    ax[n-1].set_xlabel(r"$\lambda\ (\mathring{A})$")

    if len(NameArr) == 0:
        NameArr = [""] *n
    if show_yName:
        ax[0].set_ylabel("Flux (uds)")
    ax[0].set_title(NameArr[0],loc="right", y=.5,
       rotation=270, ha="left", va="center")
    #Ploteamos los espectros
    Axe_Lined_Spectra(lambArr[0], fluxArr[0], lines, ax[0],name = NameArr[0],show_yName=True,show_xName = True)
    for i in range(1,n):   
        if show_yName:
            ax[i].set_ylabel("Flux (uds)")
        ax[i].set_title(NameArr[i],loc="right", y=.5,
           rotation=270, ha="left", va="center")
        #Ploteamos los espectros
        Axe_Lined_Spectra(lambArr[i], fluxArr[i], lines, ax[i],name = NameArr[i],show_yName=True,show_xName = False,show_LineName=False)
    return 
def Compare_Norms(defArr,normArr,fitArr = [],NameArr = False,lines = {}, title = "Spectra Normalized"):
    """
    Compare_Norms:
        Ploteamos a la izquierda como era el espectro antes de normalizarse
        con su fit, y a la derecha como lo hemos normalizado
        Podemos darlo como 
        
    Entrada:
        defArr: [normal1,normal2,...]
        con normal_i como:
            normal_i: (lambda_i,flux_1)
        normArr: [normalizado1, normalizado2,...]
        con misma estructrua que normal_i
        
        fit: [(lamb1,flux_fit1),(lamb2,flux_fit2),...]        
    """
    n = len(defArr)
    if not NameArr:
        NameArr = [False]*n
    fig, ax = plt.subplots(n,2,figsize = (60,6), sharey = False, sharex = True)
    fig.subplots_adjust(wspace = 0.1)
    fig.suptitle(title)
    if n == 1:
        lambDef,fluxDef = defArr[0]
        Axe_Blank_Spectra(lambDef, fluxDef,ax[0]) # lineas solo en el normalizado
        lambNorm,fluxNorm = normArr[0]
        Axe_Lined_Spectra(lambNorm, fluxNorm, lines,ax[1], yName = False)
        # Ponemos los fits
        if len(fitArr) != 0:
            ajustes, = ax[0].plot(fitArr[0,0],fitArr[0,1],linestyle = "dashed")# label = "ajuste")
            #ax[0].legend(handles = [ajustes], loc = "upper left")
    else:
        # En la primera ponemos las lineas
        lambDef,fluxDef = defArr[0]
        Axe_Blank_Spectra(lambDef, fluxDef, ax[0,0], show_xName=False) # Lineas solo en el normalizado
        lambNorm,fluxNorm = normArr[0]
        Axe_Lined_Spectra(lambNorm, fluxNorm,lines,ax[0,1],name = NameArr[0],show_yName = False, show_xName= False)
        
        # En la ultima ponemos los ejes
        lambDef,fluxDef = defArr[-1]
        Axe_Blank_Spectra(lambDef, fluxDef, ax[-1,0]) # Lineas solo en el normalizado
        lambNorm,fluxNorm = normArr[-1]
        Axe_Lined_Spectra(lambNorm, fluxNorm,lines,ax[-1,1],name = NameArr[-1],show_yName = False,show_LineName = False)
        for i in range(1,n-1):
            lambDef,fluxDef = defArr[i]
            Axe_Blank_Spectra(lambDef, fluxDef, ax[i,0], show_xName = False) # Lineas solo en el normalizado
            lambNorm,fluxNorm = normArr[i]
            Axe_Lined_Spectra(lambNorm, fluxNorm,lines,ax[i,1],name = NameArr[i],show_yName = False,show_xName = False, show_LineName = False)
        # Ponemos los fits
        for i in range(len(fitArr)):  
            lambFit, fluxFit = fitArr[i]
            ajustes, = ax[i,0].plot(lambFit,fluxFit,linestyle = "dashed")# label = "ajuste")
            #ax[i,0].legend(handles = [ajustes], loc = "upper right")
    
    fig.legend(ncol = ncol)
    fig.show()
    return 