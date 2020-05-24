import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import scipy.constants as const

U,N=np.genfromtxt("input/Kennlinie.dat",unpack=True)
N=unp.uarray(N,np.sqrt(N))


#Charakteristik bestimmen
def line(x,m,n):
    return m*x+n

param, cov = curve_fit(line, U[6:27] , noms(N[6:27]))
error=np.sqrt(np.diag(cov))
unparams=unp.uarray(param,error)
steigung=unparams[0]
abschnitt=unparams[1]
U_line=np.linspace(380,590)


plt.plot(U_line,param[0]*U_line+param[1],"-r",label="Plateaugerade")
plt.errorbar(U,noms(N),yerr=stds(N),fmt="xb",label="Messdaten")
plt.xlabel("Spannung $\;/\;$V")
plt.ylabel("Impulse")
plt.legend(loc="best")
plt.savefig("input/p_charakteristik.pdf")
plt.close()

#Totzeit bestimmen
N1=unp.uarray(96041,np.sqrt(96041))
N2=unp.uarray(76518,np.sqrt(76518))
N12=unp.uarray(158479,np.sqrt(158479))
N=[N1,N2,N12]
T=(N1+N2-N12)/(2*N1*N2)*120

#Zahlrohrstrom
Ni=([9837,9995,10264,10151,10184,10253,10493,11547])
Ns=unp.uarray(Ni,np.sqrt(Ni))
U_z,I_z=np.genfromtxt("input/Zaehlrohrstrom.dat",unpack=True)
I_z=unp.uarray(I_z*10**(-6),0.05*10**(-6))
Z=I_z/(const.e*Ns/60)

plt.errorbar(noms(I_z*10**(6)),noms(Z),yerr=stds(Z),fmt="x",label="Freigesetzte Ladungen")
plt.xlabel("Zählrohrstrom $\;/\;10^{-6}A$")
plt.ylabel("Freigesetzte Ladung $\;/\;$e")
plt.legend(loc="best")
plt.savefig("input/Ladungen.pdf")
plt.close()


def table():
    k=0
    print("\t c) Freigesetzte Ladungen")
    while k<8:
        print(f"{I_z[k]*10**(9)} & {Ns[k]/60} & {Z[k]} \\\\")
        k=k+1

#Ausgabe
print(f"""
Geiger-Müller-Zählrohr
\t a) Charakteristik
Gerade m={steigung}, n={abschnitt}

\t b) Totzeit
T={T}


{table()}
N:{N}
""")
