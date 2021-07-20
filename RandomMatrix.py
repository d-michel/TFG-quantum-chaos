# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as p
import random as r
import pyquaternion as pq

def randomMatrix(N, ensemble, sigmaDiag, sigmaExt):
    
    if ensemble == 'Poisson':
        RMatrix = np.zeros((N, N))
        
        for i in range(0, N):
            RMatrix[i][i] = r.gauss(0, sigmaDiag)
    
    elif ensemble == 'GOE':
        RMatrix = np.zeros((N, N))
        
        for i in range(0, N):
            RMatrix[i][i] = r.gauss(0, sigmaDiag)
            
            for j in range(i + 1, N):
                RMatrix[i][j] = r.gauss(0, sigmaExt)
                RMatrix[j][i] = RMatrix[i][j]
                
    elif ensemble == 'GUE':
        RMatrix = np.zeros((N, N), dtype = np.complex)
        
        for i in range(0, N):
            RMatrix[i][i] = r.gauss(0, sigmaDiag)
            
            for j in range(i + 1, N):
                RMatrix[i][j] = complex(r.gauss(0, sigmaExt), r.gauss(0, sigmaExt))
                RMatrix[j][i] = np.conj(RMatrix[i][j])
    
    elif ensemble == 'GSE':
        RMatrix = np.zeros((N, N), dtype = pq.Quaternion)
        for i in range(0, N):
            RMatrix[i][i] = r.gauss(0, sigmaDiag)
            
            for j in range(i + 1, N):
                RMatrix[i][j] = pq.Quaternion([r.gauss(0, sigmaExt), r.gauss(0, sigmaExt), r.gauss(0, sigmaExt), r.gauss(0, sigmaExt)])
                RMatrix[j][i] = RMatrix[i][j].conjugate
    
    else:
        raise TypeError("Ensemble type must be:\n\
            ->  'Poisson' for Poissonian Ensemble\n\
            ->  'GOE' for Gaussian Orthogonal Ensemble\n\
            ->  'GUE' for Gaussian Unitary Ensemble\n\
            ->  'GSE' for Gaussian Symplectic Ensemble")
    
    return RMatrix

def eigenvalues(RMatrix, ensemble):
    
    N = len(RMatrix)
    eigvals = np.zeros(N)
    
    for i in range(0, N):
        eigvals[i] = float(np.linalg.eigvals(RMatrix)[i])
    
    if ensemble == 'GOE' or ensemble == 'Poisson':
        eigvals = sorted(eigvals)
    
    return eigvals

def eigvalHistogram(Eigvals, n):
    
    p.hist(Eigvals, density = True, bins = n, label = 'eigenvals histogram')
    p.show()

def eigvalSpacing(eigvals):
    
    N = len(eigvals)
    spacing = np.zeros(N-1)
    
    for i in range(0, N-1):
        spacing[i] = eigvals[i+1] - eigvals[i]
    
    return spacing

def eigvalRatio(spacing):
    
    N = len(spacing) + 1
    ratio = np.zeros(N-2)
    
    for i in range(0, N-2):
        ratio[i] = spacing[i+1]/spacing[i]
    
    return ratio

def ratioHistogram(ratio, n):
    
    #fig, ax = p.subplots()
    #ax.plot(x, np.sin(x), '-b', label='Sine')
    #ax.plot(x, np.cos(x), '--r', label='Cosine')
    #ax.axis('equal')
    #leg = ax.legend();
    
    #fig, ax = p.subplots()
    p.hist(ratio, density = True, bins = n, label = 'ratios histogram')
    r = np.linspace(0,120,1200)
    Pr1 = (27/8) * (r + r**2)/(1 + r + r**2)**(5/2)
    Pr2 = (96/25) * (r + r**2)/((1 + r)**2 - 4*r/5)**(5/2)
    p.plot(r, Pr1, label='P1')
    p.plot(r, Pr2, label="P2")
    p.xlim((0, 10))
    p.legend()
    p.show()

def minRatio(spacing):
    
    N = len(spacing) + 1
    ratio1 = np.zeros(N-2)
    ratio2 = np.zeros(N-2)
    ratioMin = np.zeros(N-2)
    
    for i in range(0, N-2):
        ratio1[i] = spacing[i+1]/spacing[i]
        ratio2[i] = 1/ratio1[i]
        
        if ratio1[i] <= ratio2[i]:
            ratioMin[i] = ratio1[i]
        else:
            ratioMin[i] = ratio2[i]
    
    return ratioMin

def minRatioHistogram(ratioMin, n):
    
    p.hist(ratioMin, density = True, bins = n, label = 'min ratios histogram')
    p.xlim((0, 10))
    p.show()