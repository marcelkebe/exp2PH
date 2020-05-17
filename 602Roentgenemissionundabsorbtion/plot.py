import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from scipy.signal import find_peaks
from scipy.interpolate import UnivariateSpline
import scipy.constants

a0,N0=np.genfromtxt("input/Bragg.dat",unpack=True)

peaks,properties=find_peaks(N0,height=200)

plt.plot(a0,N0,"xk",label="Messdaten")
plt.plot(a0[peaks],N0[peaks],"or",label="max. bei $\Theta_{GM}:$"+f"{a0[22]}°"+f",\n N: {N0[22]} Imp/s")
plt.xlabel("Winkel $\Theta_{GM}\;/\;$°")
plt.ylabel("Impulsrate $N $Imp/s")
plt.legend(loc="best")
plt.savefig("plots/messdaten1.pdf")
plt.close()

h=6.62607015*10**(-34)#Js
c=3*10**8#m/s
R=h*c*1.097*10**7#1/m
e=1.602*10**(-19)
alpha = scipy.constants.value("fine-structure constant")
EKabs=8980*e
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
plt.plot(lamdba1*v,N1,"-k",label="Spektrumsverlauf")
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
aS, NS = np.genfromtxt("input/Strontium.dat",unpack=True)
azi, NZi = np.genfromtxt("input/Zirkonium.dat",unpack=True)

#Literatur Werte
Z = np.array([29, 30, 31, 35, 37, 38, 40]) #(Cu, Zn, Ga, Br, Rb, Sr, Zr)
E_klit = np.array([8.98, 9.65, 10.37, 13.47, 15.20, 16.10, 17.99])*1000 #Gleiche Reihenfolge in keV
sigma_klit=np.array([3.76,3.65,3.62,3.85,3.95,4.01,4.11])
brag=np.array([20.5,18.7,17.375,13.2,11.8,11.1,9.95])

#Berechnung von Energie, Sigma, Prozents
def getEnergiefrombragg(a):
    return h*c/(2*d*sin(a/360*2*np.pi))*6.242*10**18#eV
def plotting(x,y,label):
    plt.plot(getEnergiefrombragg(x)/1000,y,"-k")
    plt.plot(getEnergiefrombragg(x)/1000,y,"xb",label=label+" Messwerte")
    plt.xlabel("Energie $E\;/\;$keV")
    plt.ylabel("Impulsrate $N$ Imp/s")
    plt.legend(loc="best")
    plt.savefig("plots/"+f"{label}.pdf")
    plt.close()
def getsigma(E,Z):
    return Z-np.sqrt((E)/(R)-((alpha**2)*(Z**4))/(4))
def percent(a,b):
    return abs(a-b)/a*100


#IST EK und sigma berechnung
E_K=getEnergiefrombragg(brag)
sigmas=getsigma(E_K*e,Z)

#Intensitäten
I_KZmin=Nz[5]
I_KZmax=Nz[10]
I_KZ = Nz[5]+(Nz[10]-Nz[5])/2

I_KGmin=NG[1]
I_KGmax=NG[6]
I_KG = NG[1]+(NG[6]-NG[1])/2

I_KBmin=NB[2]
I_KBmax=NB[7]
I_KB = NB[2]+(NB[7]-NB[2])/2

I_KRmin=NR[3]
I_KRmax=NR[9]
I_KR = NR[3]+(NR[9]-NR[3])/2

I_KZimin=NZi[0]
I_KZimax=NZi[7]
I_KZi = NZi[0]+(NZi[7]-NZi[0])/2

I_KSmin=NS[4]
I_KSmax=NS[9]
I_KS = NS[4]+(NS[9]-NS[4])/2

#Ausgabe Sigma,EnergienK und Intensitäten
print(f"""
Sigmas\tLit\tist\t\%
CU:\t{sigma_klit[0]}\t{sigmas[0]}\t{percent(sigma_klit[0],sigmas[0])}
Zn:\t{sigma_klit[1]}\t{sigmas[1]}\t{percent(sigma_klit[1],sigmas[1])}
Ga:\t{sigma_klit[2]}\t{sigmas[2]}\t{percent(sigma_klit[2],sigmas[2])}
Br:\t{sigma_klit[3]}\t{sigmas[3]}\t{percent(sigma_klit[3],sigmas[3])}
Rb:\t{sigma_klit[4]}\t{sigmas[4]}\t{percent(sigma_klit[4],sigmas[4])}
Sr:\t{sigma_klit[5]}\t{sigmas[5]}\t{percent(sigma_klit[5],sigmas[5])}
Zi:\t{sigma_klit[6]}\t{sigmas[6]}\t{percent(sigma_klit[6],sigmas[6])}

EK\tlit\tist\t\t\%
CU:{E_klit[0]}\t{E_K[0]}\t{percent(E_klit[0],E_K[0])}
ZN:{E_klit[1]}\t{E_K[1]}\t{percent(E_klit[1],E_K[1])}
Ga:{E_klit[2]}\t{E_K[2]}\t{percent(E_klit[2],E_K[2])}
Br:{E_klit[3]}\t{E_K[3]}\t{percent(E_klit[3],E_K[3])}
Rb:{E_klit[4]}\t{E_K[4]}\t{percent(E_klit[4],E_K[4])}
Sr:{E_klit[5]}\t{E_K[5]}\t{percent(E_klit[5],E_K[5])}
Zi:{E_klit[6]}\t{E_K[6]}\t{percent(E_klit[6],E_K[6])}

Intensitäten
Ik\tImin\tImax\tBRagg
Zn:{I_KZ}\t{I_KZmin}\t{I_KZmax}\t{brag[1]}
Ga:{I_KG}\t{I_KGmin}\t{I_KGmax}\t{brag[2]}
Br:{I_KB}\t{I_KBmin}\t{I_KBmax}\t{brag[3]}
Rb:{I_KR}\t{I_KRmin}\t{I_KRmax}\t{brag[4]}
Sr:{I_KS}\t{I_KSmin}\t{I_KSmax}\t{brag[5]}
Zi:{I_KZi}\t{I_KZimin}\t{I_KZimax}\t{brag[6]}

Output
Mat\tTheta\tEK\tImin\tImax\tIK\sigma
Brom & {brag[3]} & {E_K[3]} & {I_KBmin} & {I_KBmax} & {I_KB} & {sigmas[3]}
Zink & {brag[1]} & {E_K[1]} & {I_KZmin} & {I_KZmax} & {I_KZ} & {sigmas[1]}
Gallium & {brag[2]} & {E_K[2]} & {I_KGmin} & {I_KGmax} & {I_KG} & {sigmas[2]}
Rubidium & {brag[4]} & {E_K[4]} & {I_KRmin} & {I_KRmax} & {I_KR} & {sigmas[4]}
Strontium & {brag[5]} & {E_K[5]} & {I_KSmin} & {I_KSmax} & {I_KS} & {sigmas[5]}
Zirkonium & {brag[6]} & {E_K[6]} & {I_KZimin} & {I_KZimax} & {I_KZi} & {sigmas[6]}
""")

#Plotten der verschieden MAterialien
plotting(az,Nz,"Zink")
plotting(aG,NG,"Gallium")
plotting(aB,NB,"Brom")
plotting(aR,NR,"Rubidium")
plotting(aS,NS,"Strontium")
plotting(azi,NZi,"Zirkonium")

#Sigma berenung Kupfer
n=1
m=2
l=3
EKa=8050*e#ev
EKb=8920*e#ev

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

param, covm = curve_fit(linfunc, Z[1:6], np.sqrt(E_K[1:6]/f))
error=np.sqrt(np.diag(covm))
unparams=unp.uarray(param,error)

x = np.linspace(29.9, 38.1)
plt.plot(x, param[1]*x+param[0], "--k", label=r"Lineare Regression")
plt.plot(Z[1:6], np.sqrt(E_K[1:6]/f), "xb", label=r"Messwerte")
plt.xlabel(r"Ordnungszahl $Z$")
plt.ylabel(r"$\sqrt{E_K}$ in $\sqrt{eV}$")
plt.legend(loc = "best")
plt.savefig("Moseley.pdf")
plt.close()
Ryd_lin=unparams[0]**2/h
print(f"""
Ry-Frequenz:{Ryd_lin}
Ry-Energie:{Ryd_lin*h/e}eV
Ry-Konstante:{Ryd_lin/c}""")

print(f"""
BraggPrüfung max
Winkel=28.2°\tWelll={(2*d*sin(28.2/360*2*np.pi))}\tEnergie={getEnergiefrombragg(28.2)/1000}
""")