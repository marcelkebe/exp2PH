import numpy as np
import uncertainties.unumpy as unp 
from uncertainties.unumpy import (nominal_values as noms, std_devs as stds)
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pylab import *

def funktion(x,m,n):
    return m*x+n

linie, spannung = np.genfromtxt("data.csv",delimiter=",",unpack=True)
d=(linie-1)*6

params, cov = curve_fit(funktion,d,spannung)
errors = np.sqrt(np.diag(cov))
unparams = unp.uarray(params,errors)
print(f"""
    Steigung m: {unparams[0]}
    y-Achsenabschnitt n: {unparams[1]}
    d: {d}
    """)

plt.plot(d,spannung,"xk",label="Messwerte")
plt.plot(d,funktion(d,params[0],params[1]),"-r",label="Lineare Regression")
plt.xlabel("Abstand $D\;/\;$mm")
plt.ylabel("Spannung $U\;/\;$V")
plt.legend()
plt.savefig("plot.pdf")
plt.show()