# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 09:57:12 2020

@author: Heber
"""
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
#%% valor exacto d ela derivada
up = np.cos(1.0)

h = 0.1
up_aprox = (np.sin(1+h)-np.sin(1))/h
error = up - up_aprox

print ("Valor aproximado: ",up_aprox)
print ("Valor del error: ",error)
#%%-----------------------------
# muestra 

list = [0.1, 0.01, 0.001, 0.0001, 0.00001]

aprox_values = []
errores_values = []


# aproximacion a la segunda derivada 

errores_values2 = []
aprox_values2 = []
for h in list:
    aux = (np.sin(1+h) - np.sin(1))/h
    aprox_values.append(aux)
    errores_values.append(up - aux)
    # print(h, up_aprox,error)
    # formula de segundo orden
    aux_2 = (np.sin(1+h)-np.sin(1-h))/(2*h)
    aprox_values2.append(aux_2)
    errores_values2.append(up - aux_2)

plt.loglog(list,errores_values,'o-',list,errores_values2,'o-')
plt.grid(True)
#%%---------------------
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

list = [0.1, 0.01, 0.001, 0.0001]
valor_exacto = 6*np.exp(1.0)
valor_aprox = []
valor_error = []
for h in list:
    aux = (np.exp((1+h)**2)-2*np.exp(1.0) + np.exp((1-h)**2))/h**2
    valor_aprox.append(aux)
    aux2 = abs(valor_exacto - aux)
    valor_error.append(aux2)  
plt.grid(True)   
plt.loglog(list,valor_error,'o-')

#list,valor_aprox, 'o-'


