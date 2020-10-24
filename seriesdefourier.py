# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 10:40:43 2020

@author: Código desarrollado en clase de EDP (USA - Maestría en Matemáticas Aplicadas)
"""
import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-np.pi,np.pi,500)
plt.grid()
k2 = (2/np.pi)-(4/np.pi)*((1/3)*np.cos(2*x))
k3 = (2/np.pi)-(4/np.pi)*(((1/3)*np.cos(2*x))+(1/15)*np.cos(4*x))
k4 = (2/np.pi)-(4/np.pi)*(((1/3)*np.cos(2*x))+(1/15)*np.cos(4*x))+(1/35)*np.cos(6*x)
k5 = (2/np.pi)-(4/np.pi)*(((1/3)*np.cos(2*x))+(1/15)*np.cos(4*x))+(1/35)*np.cos(6*x)+(1/63)*np.cos(8*x)
k6 = (2/np.pi)-(4/np.pi)*(((1/3)*np.cos(2*x))+(1/15)*np.cos(4*x))+(1/35)*np.cos(6*x)+(1/63)*np.cos(8*x)+(1/99)*np.cos(10*x)
plt.plot(x,k2,"-b",label="k=2")
plt.plot(x,k3,"-r",label="k=3")
plt.plot(x,k4,"r--",label="k=4")
plt.plot(x,k5,"b--",label="k=5")
plt.plot(x,k6,"g--",label="k=6")
plt.legend(loc="below left")
#-------------------------------
import numpy as np
import matplotlib.pyplot as plt
x = np.linspace(-np.pi,np.pi,100)
a0 = (16*np.pi**5)/15
d = (np.pi**2)*(k**2)-3
e = 3*np.pi*k*np.cos(np.pi*k)
c = 16/np.pi

f = np.cos(k*x)

for k in range(1,100):
    y = (a0/2)+c*(((d+e)*f)/k**5)
    
plt.plot(x,2*y)

