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
#%% Funciones
def Blank_Spectra(lamb,flux, title = "Espectro", onlyObject = False):
    """
    Blank_Spectra:
        
        Ploteamos el espectro con la menor cantidad de cambios posible.
        Perfecto para un primer acercamiento
    """
    fig, ax = plt.subplots(figsize = (15,6))
    fig.suptitle(title)
    ax.set_xlabel(r"$\lambda \ (\mathring{A})$")
    ax.set_ylabel("Flujo (unidades arbitrarias)")
    ax.plot(lamb,flux, linewidth = lWidth)
    if onlyObject:
        return fig,ax # Devuelvo los ejes sin plottear
    fig.show()
    return

def Lined_Spectra(lamb,flux,lines, title = "Lined_Spectra"):
    """
    Lined_Spectra:
        Ploteamos el espectro con las lineas que queramos marcar.
        Lines debe enviarse como un diccionario con {"line_name" : lamb,}
    """
    # Ploteamos espectro
    fig,ax = Blank_Spectra(lamb, flux, title = title, onlyObject = True)
    minLine = np.min(flux)*yscale
    maxLine = np.max(flux)*yscale
    # Ploteamos las lineas
    for name in lines:
        ax.vlines(lines[name], minLine,maxLine, label = name)
    ax.legend(loc = "best")
    fig.show()
    return

def Compare_Spectra(lambArr,fluxArr,TArr = [],lines = {}, title = "Comparison between Spectras"):
    """
    Compare_Spectra:
        Ploteamos uno al lado del otro varios espectros para compararlos. Si se dan lineas
        como parametro opcional, las coloca tambien en ambos
        
    Entrada:
        LambsArr: Array de array de longitudes de onda
        FluxArr: Array de array de Flujos
        TArr: Array de temperaturas de cada una
    """
    n = len(lambArr)
    fig,ax = plt.subplots(n,1, figsize = (60,6),sharex=True, sharey=False)
    ax[n-1].set_xlabel(r"$\lambda\ (\mathring{A})$")
    
    minLine = np.min(fluxArr)*yscale
    maxLine = np.max(fluxArr)*yscale
    
    for i in range(n):          
        ax[i].set_ylabel("Flux (uds)")
        ax[i].set_title(f"Spectra {i+1}  (Temp Estimada: {TArr[i]} K")
        #Ploteamos los espectros
        ax[i].plot(lambArr[i],fluxArr[i],linewidth = lWidth)
    
        ax[i].set_ylim(minLine,maxLine)
    for name in lines:
        ax[0].plot((lines[name],lines[name]), (minLine,maxLine),label = name) # Este para el label
        for i in range(1,n):
            ax[i].plot((lines[name],lines[name]), (minLine,maxLine))
    fig.legend()
    fig.show()
    return
    