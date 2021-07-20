# -*- coding: utf-8 -*-

import RandomMatrix as rm
import numpy as np
import matplotlib.pyplot as p

N = 100
n = 10
ensemble = 'GOE'

Eigvals = np.zeros(N*n)
Spacing = np.zeros(n*(N-1))
Ratio = np.zeros(n*(N-2))
Ratiomin = np.zeros(n*(N-2))

K = 50
# sigmaExt = np.linspace(0, 0.005, K)
sigmaExt = np.geomspace(10**(-10), 1/2, K)
sigmaExt[0] = 0
meanRatio = np.zeros(K)
meanMinRatio = np.zeros(K)

for k in range(0,len(sigmaExt)):
    for i in range(0,n):
        RMatrix = rm.randomMatrix(N,ensemble, 1, sigmaExt[k])
        eigvals = rm.eigenvalues(RMatrix, ensemble)
        spacing = rm.eigvalSpacing(eigvals)
        ratio = rm.eigvalRatio(spacing)
        ratiomin = rm.minRatio(spacing)
        j = 0
        while j < N:
            Eigvals[i*N + j] = eigvals[j]
            j += 1
            
        j = 0
        while j < (N-1):
            Spacing[i*(N-1) + j] = spacing[j]
            j += 1
            
        j = 0
        while j < (N-2):
            Ratio[i*(N-2) + j] = ratio[j]
            Ratiomin[i*(N-2) + j] = ratiomin[j]
            j += 1
            
    meanRatio[k] = np.mean(Ratio)
    meanMinRatio[k] = np.mean(Ratiomin)
    print(k)

sigmaExt[0] = sigmaExt[1]/100
sigmaExtLog = np.log(sigmaExt)

p.plot(sigmaExt, meanRatio, marker='o')
p.xlabel('Ratio')
p.show()
p.plot(sigmaExt, meanMinRatio, marker='o')
p.xlabel('Ratiomin')
p.show()
p.plot(sigmaExtLog, meanRatio, marker='o')
p.xlabel('Ratio Log')
p.show()
p.plot(sigmaExtLog, meanMinRatio, marker='o')
p.xlabel('Ratiomin Log')
p.show()