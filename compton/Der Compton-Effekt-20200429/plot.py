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
lamdbae=2*d*sin(a/360*2*np.pi)/n
print(f"lamda first: {lamdbae[0]} bis {lamdbae[-1]}")
E=h*c/lamdbae*6.242*10**18#eV
print()
peaks,properties=find_peaks(R,height=1000)
print(R[peaks])

fig, ax1 = plt.subplots()
fig.subplots_adjust(top=0.9)

color="tab:red"
ax1.set_xlabel("Wellenlänge $\lambda\;/\;pm$")
ax1.set_ylabel("Impulsrate $I_0\;/\;Impuls/s$")
ax1.plot(lamdbae/10**(-12),R,"xb",label="Messwerte")
ax1.plot(lamdbae/10**(-12),R,"-k",label="Spektrumverlauf")
ax1.set_xlim(5.605892506671436*10**(-11)/10**(-12),1.7023063582915372*10**(-10)/10**(-12))
plt.legend(loc="best")
ax1.tick_params(axis="x")

ax2=ax1.twiny()
ax2.set_xlim(8,25)
ax2.set_xlabel("Braggwinkel $\Theta\;/\;$°")
#ax2.plot(a,R)
ax1.tick_params(axis="x")

def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

#par2=ax1.twiny()
#par2.spines["top"].set_position(("axes",1.2))
#make_patch_spines_invisible(par2)
#par2.spines["top"].set_visible(True)
#par2.set_xlabel("Energie")
#par2.set_xlim()
#p3 = par2.plot(a,E)
#lines =[p3]

plt.savefig("plots/welll_int.pdf")
#plt.show()
plt.close()


plt.plot(E/1000,R,"xb",label="Messwerte")
plt.plot(E/1000,R,"-k",label="Spektrumsverlauf")
plt.plot(E[peaks[0]]/1000,R[peaks[0]],"or",label=f"$K_b=${round(E[peaks[0]]/1000,2)}keV")
plt.plot(E[peaks[1]]/1000,R[peaks[1]],"ob",label=f"$K_a=${round(E[peaks[1]]/1000,2)}keV")
plt.xlim(7.5,10)
#plt.title("Charakteristische Röntgenstrahlung")
plt.xlabel("Energie $E\;/\;keV$")
plt.ylabel("Impulsrate $I_0\;/\;Impuls/s$")
plt.legend(loc="best")
plt.savefig("plots/E_int.pdf")
#plt.show()
plt.close()

#------------------------------------------------------------------

ao,Ro=np.genfromtxt("compton/ComptonOhne.txt",unpack=True)
aA,RA=np.genfromtxt("compton/ComptonAl.txt",unpack=True)
t=200#s
No=unp.uarray(Ro,np.sqrt(Ro))
NA=unp.uarray(RA,np.sqrt(RA))
tau=90*10**(-6)#s
Io=No/(1-tau*(No))
IA=NA/(1-tau*(NA))
T=IA/Io
lamdba=2*d*sin(aA/360*2*np.pi)

def funktion(L,m,n):
    return m*L+n

params, cov = curve_fit(funktion,lamdba/10**(-12),noms(T))
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)
print(f"""
    m: {unparams[0]}
    n: {unparams[1]}
""")

plt.plot(lamdba/10**(-12),funktion(lamdba/10**(-12),params[0],params[1]),"--r",label="Lineare Regression")
plt.plot(lamdba/10**(-12),noms(T),"xb",label="Messwerte")
#plt.title("Transmission $T(\lambda)$")
plt.xlabel("Wellenlänge $\lambda\;/\;$pm")
plt.ylabel("Transmission $T\;/\;$%")
plt.legend(loc="best")
#plt.show()
plt.savefig("plots/trans.pdf")
plt.close()


#--------------------------------------------------------------------
t=300#s
me=9.109*10**(-31)#kg
I0=unp.uarray(2731,np.sqrt(2731))#Impulse 
I1=unp.uarray(1180,np.sqrt(1180))#Impulse ungestreut
I2=unp.uarray(1024,np.sqrt(1024))#Impulse gestreut
T1=I1/I0
T2=I2/I0
print(f"""
    T1ung: {T1}
    T2g: {T2}
""")

def getlambda(T):
    return (T-unparams[1])/unparams[0]#pm
print(f"""
    Wellenlängen
    lambda1 (ungest. R.)\t ({getlambda(T1)})\tpm
    lambda2 (gest. R.)\t\t ({getlambda(T2)})\tpm
    delta \t\t\t ({getlambda(T2)-getlambda(T1)})\tpm
    fehler \t\t\t ({(abs(2.4-(getlambda(T2)-getlambda(T1)))/2.4)*100})
    
    I=0 {I0}
    Iung {I1}
    Ig {I2}
""")
print(f"""
    #Char. Strahlung
    Kb: {round(E[peaks[0]]/1000,2)} keV
    TKb: {round(a[peaks[0]],2)}°
    LKb: {round(lamdbae[peaks[0]]/10**(-12),2)}pm
    Ka: {round(E[peaks[1]]/1000,2)} keV
    TKa: {round(a[peaks[1]],2)}°
    LKa: {round(lamdbae[peaks[1]]/10**(-12),2)}pm
    """)

k=len(lamdba)
x=0
while x<=k-1:
    print(f"""{ao[x]} & {lamdba[x]} & {No[x]} & {noms(Io[x])} & {NA[x]} & {noms(IA[x])} & {noms(T[x])}""")
    x+=1

