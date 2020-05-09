import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline

a0,N0=np.genfromtxt("input/Bragg.dat",unpack=True)

peaks,properties=find_peaks(N0,height=200)

plt.plot(a0,N0,"xk",label="Messdaten")
plt.plot(a0[peaks],N0[peaks],"or",label="max. bei $\Theta_{GM}:$"+f"{a0[22]}°"+f",\n N: {N0[22]} Imp/s")
plt.xlabel("Winkel $\Theta_{GM}\;/\;$°")
plt.ylabel("Impulsrate $N $Imp/s")
plt.legend(loc="best")
plt.savefig("plots/messdaten1.pdf")
plt.close()

n=1
d=201.4*10**(-12)
v=10**(12)#faktor for m in pm
a1,N1=np.genfromtxt("input/Emissionsspektrum.dat",unpack=True)
lamdba1=2*d*sin(a1/360*2*np.pi)/n
peaks,properties=find_peaks(N1,height=1000)
#print(f"peaks:{peaks}\n({lamdba1[peaks]}/{N1[peaks]})")

a1_1=a1[0:130]
N1_1=N1[0:130]
lamdba1_1=lamdba1[0:130]

a1_2=a1[130:len(a1)]
N1_2=N1[130:len(a1)]
lamdba1_2=lamdba1[130:len(a1)]

#FWHM
spline=UnivariateSpline(lamdba1_1,N1_1-np.max(N1_1)/2,s=0)
r1_1, r1_2=spline.roots()*v
spline=UnivariateSpline(lamdba1_2,N1_2-np.max(N1_2)/2,s=0)
r2_1, r2_2=spline.roots()*v

#PLOTS - Wellenlänge
plt.plot(lamdba1*v,N1,"xb",label="Messwerte")
plt.plot(lamdba1*v,N1,"-k")
plt.axvspan(r1_1,r1_2,facecolor="g",alpha=0.5)
plt.axvspan(r2_1,r2_2,facecolor="y",alpha=0.5)
plt.xlabel("$\lambda\;/\;$pm")
plt.ylabel("Impulsrate $N$ Imp/s")
plt.legend(loc="best")
plt.show()
plt.savefig("plots/spektrum.pdf")

print(f"""
Halbwertsbreiten FWHM
Ka\t{round(r2_2-r1_2,2)}pm
Kb\t{round(r1_2-r1_1,2)}pm""")

#PLOTS - Energien
h=6.62607015*10**(-34)#Js
c=3*10**8#m/s    
def getEnergie(l):
    return h*c/l*6.242*10**18#eV

h=6.62607015*10**(-34)#Js
c=3*10**8#m/s
    
E=h*c/lamdba1*6.242*10**18#eV
print(f"""Roots
1:r1={r1_1} r2={r1_2}
2:r1={r2_1} r2={r2_2}""")
print(f"""
Eapha: {round(E[peaks[0]]/1000,2)} keV
Ebeta: {round(E[peaks[1]]/1000,2)} keV
""")

E_da=getEnergie((r2_2-r2_1)*1/v)
E_db=getEnergie((r1_2-r1_1)*1/v)
print(f"""
E_da:{E_da} eV
E_db:{E_db} eV""")
A_a=E[peaks[0]]/E_da
A_b=E[peaks[1]]/E_db
print(f"""
Auflösungsvermögen A
A_a: {A_a}
A_b: {A_b}""")

az, Nz = np.genfromtxt("input/Zink.dat",unpack=True)
aG, NG = np.genfromtxt("input/Gallium.dat",unpack=True)
aB, NB = np.genfromtxt("input/Brom.dat",unpack=True)
aR, NR = np.genfromtxt("input/Rubidium.dat",unpack=True)

def plotting(x,y,label):
    plt.plot(x,y,"-k")
    plt.plot(x,y,"xb",label=label+" Messwerte")
    plt.xlabel("Winel $\Theta\;/\;$°")
    plt.ylabel("Impulsrate $N$ Imp/s")
    plt.legend(loc="best")
    plt.savefig("plots/"+f"{label}.pdf")
    #plt.show()

plotting(az,Nz,"Zink")
plotting(aG,NG,"Gallium")
plotting(aB,NB,"Brom")
plotting(aR,NR,"Rubidium")