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
#plt.show()
plt.savefig("plots/spektrum.pdf")
plt.close()

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
print(f"""Roots (lambda/pm)
1:r1={r1_1} r2={r1_2}
2:r1={r2_1} r2={r2_2}""")
print(f"""
Ebeta: {round(E[peaks[0]]/1000,2)} keV
Ealpha: {round(E[peaks[1]]/1000,2)} keV
""")

E_da=getEnergie(r2_1*1/v)-getEnergie(r2_2*1/v)
E_db=getEnergie(r1_1*1/v)-getEnergie(r1_2*1/v)
print(f"ROOT(Energie/eV):\npeakb: {getEnergie(r1_1*1/v)},{getEnergie(r1_2*1/v)},\npeaka: {getEnergie(r2_1*1/v)},{getEnergie(r2_2*1/v)}")
print(f"""
E_da:{E_da} eV
E_db:{E_db} eV""")
A_a=E[peaks[1]]/E_da
A_b=E[peaks[0]]/E_db
print(f"""
Auflösungsvermögen A
A_a: {A_a}
A_b: {A_b}""")

az, Nz = np.genfromtxt("input/Zink.dat",unpack=True)
aG, NG = np.genfromtxt("input/Gallium.dat",unpack=True)
aB, NB = np.genfromtxt("input/Brom.dat",unpack=True)
aR, NR = np.genfromtxt("input/Rubidium.dat",unpack=True)

def getEnergiefrombragg(a):
    return h*c/(2*d*sin(a/360*2*np.pi))*6.242*10**18#eV

def plotting(x,y,label):
    plt.plot(getEnergiefrombragg(x)/1000,y,"-k")
    plt.plot(getEnergiefrombragg(x)/1000,y,"xb",label=label+" Messwerte")
    plt.xlabel("Energie $E\;/\;$keV")
    plt.ylabel("Impulsrate $N$ Imp/s")
    plt.legend(loc="best")
    plt.savefig("plots/"+f"{label}.pdf")
    plt.show()

plotting(az,Nz,"Zink")
plotting(aG,NG,"Gallium")
plotting(aB,NB,"Brom")
plotting(aR,NR,"Rubidium")

e=1.602*10**(-19)
EKabs=8980*e
n=1
m=2
l=3
EKa=8050*e#ev
EKb=8920*e#ev
R=1.097*10**7#1/m

Z = np.array([29, 30, 31, 35, 37, 38, 40]) #(Cu, Zn, Ga, Br, Rb, Sr, Zr)
E_k = np.array([8.98, 9.65, 10.37, 13.47, 15.20, 16.10, 17.99]) *1000 #Gleiche Reihenfolge in keV
#sigma1=Z[0]+np.sqrt(EKabs/R)
#sigma2=1/(n**2*R)*(n**2*R*Z[0]+np.sqrt(m**2*n**2*R*(EKa*n**2-R*sigma1**2+2*R*sigma1**2*Z[0]-R*Z[0]**2)))
#sigma3=1/(n**2*R)*(n**2*R*Z[0]+np.sqrt(l**2*n**2*R*(EKb*n**2-R*sigma1**2+2*R*sigma1**2*Z[0]-R*Z[0]**2)))

sigma1 = Z[0]-np.sqrt((EKabs)/(R))
sigma2 = Z[0]-np.sqrt( 4*(Z[0]-sigma1)**2-(EKa)/(R)*4 )
sigma3 = Z[0]-np.sqrt( 9*(Z[0]-sigma1)**2-(EKb)/(R)*9 )
print(f"""
sigma1={sigma1}
sigma2={sigma2}
sigma3={sigma3}""")


E_kRb = getEnergiefrombragg(11.76)
E_kSr=getEnergiefrombragg(11.09)
E_kZn=getEnergiefrombragg(18.65)
E_kGa=getEnergiefrombragg(17.34)
E_kBr=getEnergiefrombragg(13.2)
E_ks = np.array([E_kZn, E_kGa, E_kBr, E_kRb, E_kSr])
f=6.242*10**18
def linfunc(a,b,x):
    return a*x+b

param, covm = curve_fit(linfunc, Z[1:6], np.sqrt(E_ks/f))
x = np.linspace(29.9, 38.1)
plt.plot(x, param[1]*x+param[0], "--k", label=r"Lineare Regression")
plt.plot(Z[1:6], np.sqrt(E_ks/f), "xb", label=r"Messwerte")
plt.xlabel(r"Ordnungszahl $Z$")
plt.ylabel(r"$\sqrt{E_K}$ in $\sqrt{eV}$")
plt.legend(loc = "best")
plt.savefig("Moseley.pdf")
plt.show()