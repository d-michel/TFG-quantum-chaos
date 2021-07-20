# -*- coding: utf-8 -*-

import RandomMatrix as rm
import numpy as np

N = 10
n = 10
ensemble = 'GOE'

Eigvals = np.zeros(N*n)
Spacing = np.zeros(n*(N-1))
Ratio = np.zeros(n*(N-2))
Ratiomin = np.zeros(n*(N-2))

for i in range(0,n):
    RMatrix = rm.randomMatrix(N,ensemble, 1, 1/2)
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

rm.eigvalHistogram(Eigvals, 2*n)
rm.ratioHistogram(Ratio, 2*N**2)
rm.minRatioHistogram(Ratiomin, N)
meanRatio = np.mean(Ratio)
meanMinRatio = np.mean(Ratiomin)
print(meanRatio, meanMinRatio)