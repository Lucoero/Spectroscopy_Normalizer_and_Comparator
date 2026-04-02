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
ncol = 1 # Numero de columnas para la leyenda
#%% Funciones
def Pad_Array(arr): # Gracias, https://stackoverflow.com/questions/24494356/how-to-find-min-max-values-in-array-of-variable-length-arrays-with-numpy   
    nArr = len(arr)
    M = max(len(a) for a in arr) # Calculamos la longitud maxima de los arrays
    out = np.zeros((nArr,M),dtype = object)
    for i in range(nArr):
        out[i] = np.array(list(arr[i])+[np.nan]*(M-len(arr[i])))
    return out


def Blank_Spectra(lamb,flux, title = "Espectro"):
    """
    Blank_Spectra:
        
        Ploteamos el espectro con la menor cantidad de cambios posible.
        Perfecto para un primer acercamiento
    """
    fig, ax = plt.subplots(figsize = (15,6))
    fig.suptitle(title)
    Axe_Blank_Spectra(lamb, flux, ax)
    fig.show()
    return

def Axe_Blank_Spectra(lamb,flux,ax):
    ax.set_xlabel(r"$\lambda \ (\mathring{A})$")
    ax.set_ylabel("Flujo (unidades arbitrarias)")
    ax.plot(lamb,flux, linewidth = lWidth)
    return 
def Lined_Spectra(lamb,flux,lines, title = "Lined_Spectra"):
    """
    Lined_Spectra:
        Ploteamos el espectro con las lineas que queramos marcar.
        Lines debe enviarse como un diccionario con {"line_name" : lamb,}
    """
    # Ploteamos espectro
    fig,ax = plt.subplots(figsize = (15,6))
    fig.suptitle(title)
    Axe_Lined_Spectra(lamb,flux,lines,ax)
    fig.legend(ncol = ncol)
    fig.show()
    return
 
def Axe_Lined_Spectra(lamb,flux,lines,ax):
    Axe_Blank_Spectra(lamb, flux, ax)
    minLine = np.min(flux)*yscale
    maxLine = np.max(flux)*yscale
    # Ploteamos las lineas
    for name in lines:
        ax.plot((lines[name],lines[name]), (minLine,maxLine),label = name,linestyle = "dashed",linewidth = lWidth * lineScale) # Este para el label
        ax.annotate(name,xy = (lines[name] + 5,maxLine -0.15*maxLine), xycoords = "data", ha = 'center', va = 'bottom',rotation = 'vertical')
    return 
    
def Compare_Spectra(lambArr,fluxArr,TArr = [],lines = {}, title = "Comparison between Spectras",show_T = True):
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
    Axe_Compare_Spectra(lambArr,fluxArr,ax, TArr = TArr, lines = lines,show_T = show_T)
    fig.legend(ncol = ncol)
    fig.show()
    return

def Axe_Compare_Spectra(lambArr,fluxArr, ax,TArr = [],lines = {},show_T = True):
    """
    Axe_Compare_Spectra:
        Hace lo que hace Compare_Spectra, pero para utilizarse por otros plots
    """
    n = len(ax)
    ax[n-1].set_xlabel(r"$\lambda\ (\mathring{A})$")
    
    # Esto da problemas si el array de flujos no es una matriz cuadrada
        # Arreglo: creamos una matriz cuadrada con nans y buscamos ahi
    pad = Pad_Array(fluxArr) 
    minLine = max(np.nanmin(pad)*yscale,0)
    maxLine = min(np.nanmax(pad)*yscale,5)
    if len(TArr) == 0:
        TArr = ["Not Estimated"] *n
    for i in range(n):          
        ax[i].set_ylabel("Flux (uds)")
        if show_T:
            ax[i].set_title(f"{i+1}  (T: {TArr[i]} K)",loc="right", y=.5,
               rotation=270, ha="left", va="center")
        #Ploteamos los espectros
        ax[i].plot(lambArr[i],fluxArr[i],linewidth = lWidth)
    
        ax[i].set_ylim(minLine,maxLine)
    for name in lines:
        ax[0].plot((lines[name],lines[name]), (minLine,maxLine),label = name,linestyle = "dashed",linewidth = lWidth * lineScale) # Este para el label
        ax[0].annotate(name,xy = (lines[name] + 5,maxLine -0.15*maxLine), xycoords = "data", ha = 'center', va = 'bottom',rotation = 'vertical')
        for i in range(1,n):
            ax[i].plot((lines[name],lines[name]), (minLine,maxLine), linestyle = "dashed",linewidth = lWidth * lineScale)
    return 
def Compare_Norms(defArr,normArr,fitArr = [],TArr = [],lines = {}, title = "Spectra Normalized", onlyObject = False):
    """
    Compare_Norms:
        Ploteamos a la izquierda como era el espectro antes de normalizarse
        con su fit, y a la derecha como lo hemos normalizado
        Podemos darlo como 
        
    Entrada:
        defArr: [normal1,normal2,...]
        con normal_i como:
            normal_i: [lambda_i,flux_1]
        normArr: [normalizado1, normalizado2,...]
        con misma estructrua que normal_i
        
        fit: [[lamb1,flux_fit1],[lamb2,flux_fit2],...]]        
    """
    n = len(defArr)
    fig, ax = plt.subplots(n,2,figsize = (60,6), sharey = False, sharex = True)
    fig.subplots_adjust(wspace = 0.1)
    fig.suptitle(title)
    if n == 1:
        Axe_Blank_Spectra(defArr[0,0], defArr[0,1],ax[0]) # lineas solo en el normalizado

        Axe_Lined_Spectra(normArr[0,0], normArr[0,1], lines,ax[1])
        # Ponemos los fits
        ajustes, = ax[0].plot(fitArr[0,0],fitArr[0,1],linestyle = "dashed")# label = "ajuste")
        #ax[0].legend(handles = [ajustes], loc = "upper left")
    else:
        Axe_Compare_Spectra(defArr[:,0], defArr[:,1], ax[:,0],TArr = TArr, lines = {},show_T = False) # Lineas solo en el normalizado
        Axe_Compare_Spectra(normArr[:,0], normArr[:,1],ax[:,1],TArr = TArr, lines = lines)
        # Ponemos los fits
        for i in range(len(fitArr)):  
            ajustes, = ax[i,0].plot(fitArr[i,0],fitArr[i,1],linestyle = "dashed")# label = "ajuste")
            #ax[i,0].legend(handles = [ajustes], loc = "upper right")
    
    fig.legend(ncol = ncol)
    fig.show()
    return 