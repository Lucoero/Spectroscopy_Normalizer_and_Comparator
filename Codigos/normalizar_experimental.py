import numpy as np
from scipy.signal import find_peaks, peak_prominences, peak_widths, find_peaks_cwt, savgol_filter

import normalizar
from parametros import *
import Herramientas

def picos2(wave, flux, lineas, d):
    pki = []
    for li in lineas:
        center = np.searchsorted(wave,li,side='right')
        a = center - d
        b = center + d
        rango = -1*flux[a:b]
        ind, prop = find_peaks(rango, distance = 2*d)
        if len(ind)==0:
            print("error")
        else:
            ind +=a
            pki.append(*ind)
    return pki


def continuo(fr,pr=0.1,d=5,sg=False):
    if sg:
        f = savgol_filter(fr,5,2)
    else:
        f = fr
    pk1, _ = find_peaks(-1*f,prominence=pr,distance=d)
    w1 = peak_widths(-1*f,pk1,rel_height=0.5)
    pk2, _ = find_peaks(f,prominence=pr,distance=d)
    w2 = peak_widths(f,pk2,rel_height=0.5)
    pks = np.append(pk1,pk2)
    ws = np.append(w1,w2)
    aux = np.zeros_like(fr)
    for i in range(len(pks)):
        a = int(pks[i] - ws[i])
        b = int(pks[i] + ws[i])
        lsup = f[b]
        linf = f[a]
        pe = (lsup - linf)/(2*ws[i])
        for t in range(int(a),int(b+1)):
            aux[t] = linf + pe*(t-a)
    for j in range(len(f)):
        if aux[j] == 0:
            aux[j] = f[j]
    return aux

def Norm_Agg(l,fr,rang,pr=0.1,d=10):

    a = continuo(fr,sg=False) # Sg es usar savgol o no

    norm = fr/a

    plt.plot(l,norm)
    plt.show()
    return


def npcont(l,f,pr=0.1,w=5,d=3):
    pk1, _ = find_peaks(-1*f,prominence=pr,width=w,distance=d)
    pk2, _ = find_peaks(f,prominence=pr,width=w,distance=d)
    pks = np.append(pk1,pk2)
    plt.plot(l,f,l[pks],f[pks],'*')
    plt.show()
    return pks

def diflin(f,p,d): #aquí p es el índice del pico
    aux = np.zeros(len(f))
    for i in range(len(p)):
        a = p[i] - d
        b = p[i] + d
        lsup = f[b]
        linf = f[a]
        pe = (lsup - linf)/(2*d)
        for t in range(int(a),int(b+1)):
            aux[t] = linf + pe*(t-a)
    for j in range(len(f)):
        if aux[j] == 0:
            aux[j] = f[j]
    return aux




def hola(wave,flux,lineas, d = 15):
    pks = picos2(wave,flux,lineas,d)
    rr = diflin(flux,pks,d)
    norm = flux/rr
    return norm
