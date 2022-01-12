''' PC2 - Conduction dans un anneau cylindrique - Résolution numérique
    ENSTA Paris - MF202
    Élève : Pedro Morel Rosa
    06/01/2021 '''
    
# Importation des modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tck
import seaborn as sns
from scipy import fft, ifft
sns.set()

# Définition des paramètres
a = 3
b = 4
M = 100
N = 160

# Conditions aux limites
'''Les deux conditions à la limite sont définies ci-dessous.
Pour obtenir l'autre image, décommentez simplement 
les deux derniers Tb et Ta.'''
Ta = 40 + 30*np.cos(np.arange(0,1,1/N)*np.pi*8)
Tb = 90 + 10*np.sin(np.arange(0,1,1/N)*np.pi*2)
#Tb = [80]*N
#Ta = [30]*N

# Valeurs de r et theta discretisés
r_j = list(a + (np.arange(1, M+1) - 1) * (b - a)/(M - 1))
theta_l = 2 * np.pi * (np.arange(1, N+1) - 1)/N

# Listes crées pour organiser les modes de l'espace spectral
nfft = np.concatenate((np.arange(0, N/2+1),np.arange(-N/2+1, 0)))
nphys = np.arange(0,N)

# Résolutin numérique
Tn = []
dr = (b-a)/(M-1)

for n in nphys:
    A = np.zeros((M,M))
    B = np.zeros(M)
    
    A[0,0] = 1
    A[M-1,M-1] = 1
    B[0] = fft(Ta)[n]
    B[-1] = fft(Tb)[n]
    for i in range(1, M-1):
        A[i, i-1] = - 1/(2*r_j[i+1]*dr) + 1/(dr**2)
        A[i, i] = -2/(dr**2) - nfft[n]**2/(r_j[i+1]**2)
        A[i, i+1] = 1/(2*r_j[i+1]*dr) + 1/(dr**2)
    
    Tn.append(np.linalg.solve(A,B))

Tn = list(map(list,Tn))


# Application de la transformé inverse
Tn3 = []
for i in range(M):
    Tn2 = []
    for j in range(N):
        Tn2.append(float(Tn[j][i]))
    Tn3.append(list(ifft(Tn2)))
    
# Organisation des résultats
TN = []
for i in range(N):
    Tn2 = []
    for j in range(M):
        Tn2.append(float(Tn3[j][i]))
    TN.append(Tn2)
    
# Obtention de l'image
space_theta = np.radians(np.linspace(0, 360, N))
space_r = np.linspace(2, 3, M)

r, theta = np.meshgrid(space_r, space_theta)
theta = theta/np.pi

sns.set()
fig, axes = plt.subplots(dpi = 140)
plt.contourf(theta, r , TN, 200, cmap = 'rainbow')
cbar = plt.colorbar()
cbar.set_label('Temperature', rotation=90)
axes.xaxis.set_major_formatter(tck.FormatStrFormatter('%g $\pi$'))
axes.xaxis.set_major_locator(tck.MultipleLocator(base=0.5))
plt.xlabel("Contour du Cylindre")
plt.ylabel("Distance Radiale")
plt.show()
