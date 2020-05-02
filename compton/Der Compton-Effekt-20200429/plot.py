import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *
from scipy.signal import find_peaks


d=201.4*10**(-12)
n=1
h=6.62607015*10**(-34)#Js
c=3*10**8#m/s
a, R = np.genfromtxt("EmissionCu.dat",unpack=True)
lamdba=2*d*sin(a/360*2*np.pi)/n

E=h*c/lamdba*6.242*10**18#eV

peaks,properties=find_peaks(R,height=1000)
print(R[peaks])

plt.plot(lamdba/10**(-12),R)
plt.title("Spektrum der Kupfer-Röntgenröhre")
plt.xlabel("Wellenlänge $\lambda\;/\;pm$")
plt.ylabel("Intensität $I_0\;/\;Impuls/s$")
plt.savefig("plots/welll_int.pdf")
#plt.show()
plt.close()


plt.plot(E/1000,R,"k")
plt.plot(E[peaks[0]]/1000,R[peaks[0]],"xr",label=f"$K_b=${round(E[peaks[0]]/1000,2)}keV")
plt.plot(E[peaks[1]]/1000,R[peaks[1]],"xb",label=f"$K_a=${round(E[peaks[1]]/1000,2)}keV")
plt.xlim(7.5,10)
plt.title("Charakteristische Röntgenstrahlung")
plt.xlabel("Energie $E\;/\;keV$")
plt.ylabel("Intensität $I_0\;/\;Impuls/s$")
plt.legend()
plt.savefig("plots/E_int.pdf")
#plt.show()
plt.close()

#------------------------------------------------------------------

ao,Ro=np.genfromtxt("compton/ComptonOhne.txt",unpack=True)
aA,RA=np.genfromtxt("compton/ComptonAl.txt",unpack=True)
t=200#s
No=unp.uarray(Ro*t,np.sqrt(Ro*t))
NA=unp.uarray(RA*t,np.sqrt(RA*t))
tau=90*10**(-9)#s
Io=No/(1-tau*(No))
IA=NA/(1-tau*(NA))
T=IA/Io
lamdba=2*d*sin(aA/360*2*np.pi)

def funktion(L,m,n):
    return m*L+n

params, cov = curve_fit(funktion,lamdba/10**(-12),noms(T))
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)

plt.plot(lamdba/10**(-12),funktion(lamdba/10**(-12),params[0],params[1]),"--r",label="Lineare Regression")
plt.plot(lamdba/10**(-12),noms(T),"b",label="Messwerte")
plt.title("Transmission $T(\lambda)$")
plt.xlabel("Wellenlänge $\lambda\;/\;$pm")
plt.ylabel("Transmission $T\;/\;$%")
plt.legend(loc="best")
#plt.show()
plt.savefig("plots/trans.pdf")
plt.close()


#--------------------------------------------------------------------
t=300#s
me=9.109*10**(-31)#kg
I0=2731#Impulse
I1=1180#Impulse
I2=1024#Impulse
T1=I1/I0
T2=I2/I0

def getlambda(T):
    return T-unparams[1]/unparams[0]#pm
print(f"""
    Wellenlängen
    lambda1 (ungest. R.)\t ({getlambda(T1)})pm
    lambda2 (gest. R.)\t\t ({getlambda(T2)})pm
    delta \t\t ({getlambda(T2)-getlambda(T1)})pm
    delta vergl. \t\t ({round(h/me*c,2)})pm
""")

