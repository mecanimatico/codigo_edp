# -*- coding: utf-8 -*-
"""
Created on May 2020

@authors: Rafael Torres Pabon & Heber Pareja Reinoso

University: Sergio Arboleda (U.S.A.)
Bogota Colombia
"""
#%% Se importan las librerias y paquetes necesarios
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
#%% Se establecen los parametros del modelo con sus valores
beta = 0.514
betaE = 0.25
betaI = 1
betaV=0.9
sigma=0.5
sigma_menos_1 = 2
gamma=0.2
gamma_menos_1 = 5
delta_menos_1 = 365
delta = 1/365
mu = 0.000000055
r = 0.0000714
k = 0.0001857
alfa = 0.0000093
teta_menos_1= 365
phi=1/20 
theta = 1/365
d1 = 0.05
d2 = 0.05
d3 = 0.025
d4 = 0.001
d5 = 0.0
N = 60
l = 6
Delta_x = l/N
Delta_t1 = 0.05
Delta_t2 = 0.05
Delta_t3 = 0.15
Delta_t4 = 4.5
Delta_t5 = 0
Aux_Delta_t=np.array([Delta_t1,Delta_t2,Delta_t3,Delta_t4])
delta_t=min(Aux_Delta_t)
#%% Se define la condicion inicial en el eje x
x = np.linspace(-3,3,N+1)
def aux_S3(x):
    if -3 < x or -3==x:
        return 0.86*(np.exp(-(x/(1.4))**2))
def aux_E3(x):
    if -3 < x or -3==x:
        return 0
def aux_I3(x):
    if -3 < x or -3==x:
        return 0.04*(np.exp(-(x/(1.4))**2))
def aux_R3(x):
    if -3 < x or -3==x:
        return 0
def aux_V3(x):
    if -3 < x or -3==x:
        return 0.1*(np.exp(-(x/(1.4))**2))
#%% Se grafica la condicion inicial
fig = plt.figure(facecolor='w')
ax = fig.add_subplot(111, facecolor='#dce3de', axisbelow=True)
ax.plot(x,np.vectorize(aux_S3)(x), 'b', alpha=0.5, lw=2, label='S')
ax.plot(x,np.vectorize(aux_E3)(x), 'cyan', alpha=0.5, lw=2, label='E')
ax.plot(x,np.vectorize(aux_I3)(x), 'r', alpha=0.5, lw=2, label='I')
ax.plot(x,np.vectorize(aux_R3)(x), 'g', alpha=0.5, lw=2, label='R')
ax.plot(x,np.vectorize(aux_V3)(x), '-', alpha=0.5, lw=2, label='V')
plt.title("Condición inicial SEIRV para la Influenza")
ax.set_xlabel('Espacio')
ax.set_ylabel('Tasa poblacional')
ax.set_ylim(-0.01,1.2)
ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', lw=1, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
#%%Se calcula la matriz como una funcion para cada ecuacion con difusion
def calcularmatriz (N, Delta_x, Delta_t, d):
    miu = (d * Delta_t)/(Delta_x ** 2)
    A = np.zeros([N+1,N+1])
    np.fill_diagonal(A,1 - 2*miu)
    
    for i in range(1,N+1):
        A[i,i-1] = miu
        A[i-1,i] = miu
    
    A[0][1] =2 *miu
    A[N][N-1]=2*miu
    return A
#%%# Se calcula la matriz de difusion para cada ecuacion
# Matriz para los suceptibles, a S le corresponde d1= 0.05
A1 = calcularmatriz(N, Delta_x, delta_t, d1)
# Matriz para los vacunados, a V le corresponde d2 = 0.05
A2 = calcularmatriz(N, Delta_x, delta_t, d2)
# Matriz para los expuestos, a E le corresponde d3 = 0.025
A3 = calcularmatriz(N, Delta_x, delta_t, d3)
# Matriz para los infectados, a I le corresponde d4 = 0.001
A4 = calcularmatriz(N, Delta_x, delta_t, d4)
# Matriz para los infectados, a R le corresponde d5 = 0.0
A5 = calcularmatriz(N, Delta_x, delta_t, d5)
#%% Se crean los coeficientes de reaccion para cada una de las ecuaciones
def f( S, E, I, R, V):
    fS = (-beta*betaE*E*S)-(beta*betaI*I*S)+(alfa*I*S)-(phi*S)-(r*S)
    return fS

def g( S, E, I, R, V):
     gE = (beta*betaE*E*S)+(beta*betaI*I*S)+(beta*betaE*betaV*E*V)
     +(beta*betaI*betaV*I*V)
     return gE

def h( S, E, I, R, V):
     hI =  (sigma*E)-(r+alfa+gamma)*I+(alfa*(I**(2)))
     return hI

def m( S, E, I, R, V):
     mR =  (k*E)+(gamma*I)-(r*R)-(delta*R)+(alfa*I*R)
     return mR
 
def n( S, E, I, R, V):
     nV = (-beta*betaE*betaV*E*V)-(beta*betaI*betaV*I*V)+(alfa*I*V)-(r*V)
     return nV
#%%Se establece el tiempo final, teniendo en cuenta el numero de dias
tiempofinal = 60
nt=int(tiempofinal/delta_t)
#%% Se calcula por diferencias centradas cada ecuacion del modelo
S_new=[np.vectorize(aux_S3)(x)]
V_new=[np.vectorize(aux_V3)(x)]
E_new=[np.vectorize(aux_E3)(x)]
I_new=[np.vectorize(aux_I3)(x)]
R_new=[np.vectorize(aux_R3)(x)]
for i in range(nt):
    S_new.append(S_new[i]+
                 delta_t*f( S_new[i], E_new[i],I_new[i],
                            R_new[i], V_new[i]))

    V_new.append(V_new[i]+delta_t*n(S_new[i], E_new[i],I_new[i],
                            R_new[i], V_new[i]))

    E_new.append(E_new[i]+delta_t*g( S_new[i], E_new[i],I_new[i],
                            R_new[i], V_new[i]))

    I_new.append(A4.dot(I_new[i])+
                 delta_t*h( S_new[i], E_new[i],I_new[i],
                            R_new[i], V_new[i]))

    R_new.append(A5.dot(R_new[i])+delta_t*m( S_new[i], E_new[i],I_new[i],
                            R_new[i], V_new[i]))    
#%% Se realizan las graficas de las soluciones de SVEIR
fig=plt.figure()
ax=plt.gca()
nplot=11
ns=int(6/delta_t)
#%%Se grafica la solucion de la poblacion de Suceptibles
for i in range(nplot):
    plt.grid(True)
    plt.title("S(x,t) versus x, t=0,6,...,60.")
    plt.xlabel('x')
    plt.ylabel("S(x,t) t=0,6,...,60.")
    plt.plot(x,S_new[i*ns])
## Escribir en la terminal %matplotlib auto
fig=plt.figure()
ax=plt.gca()
def actualizar(i):
    ax.clear()
    plt.plot(x,S_new[i*10],"blue")
    plt.title("S(x,t) versus x, t=0,6,...,60.")
    plt.xlim(-3,3)
    plt.ylim(0,1)
ani=animation.FuncAnimation(fig,actualizar,range(int(len(S_new)/10)))
plt.show() 
#%%Se grafica la solucion de la poblacion de Vacunados
for i in range(nplot):
    plt.grid(True)
    plt.title("V(x,t) versus x, t=0,6,...,60.")
    plt.xlabel('x')
    plt.ylabel("V(x,t) t=0,6,...,60.")
    plt.plot(x,V_new[i*ns])
## Escribir en la terminal %matplotlib auto
fig=plt.figure()
ax=plt.gca()
def actualizar(i):
    ax.clear()
    plt.plot(x,V_new[i*10],"blue")
    plt.title("V(x,t) versus x, t=0,6,...,60.")
    plt.xlim(-3,3)
    plt.ylim(0,0.11)
ani=animation.FuncAnimation(fig,actualizar,range(int(len(V_new)/10)))
plt.show() 
#%%Se grafica la solucion de la poblacion de Expuestos
for i in range(nplot):
    plt.grid(True)
    plt.title("E(x,t) versus x, t=0,6,...,60.")
    plt.xlabel('x')
    plt.ylabel("E(x,t) t=0,6,...,60.")
    plt.plot(x,E_new[i*ns])
## Escribir en la terminal %matplotlib auto
fig=plt.figure()
ax=plt.gca()
def actualizar(i):
    ax.clear()
    plt.plot(x,E_new[i*10],"blue")
    plt.title("E(x,t) versus x, t=0,6,...,60.")
    plt.xlim(-3,3)
    plt.ylim(0,0.6)
ani=animation.FuncAnimation(fig,actualizar,range(int(len(E_new)/10)))
plt.show()
#%%Se grafica la solucion de la poblacion de Infectados
for i in range(nplot):
    plt.grid(True)
    plt.title("I(x,t) versus x, t=0,6,...,60.")
    plt.xlabel('x')
    plt.ylabel("I(x,t) t=0,6,...,60.")
    plt.plot(x,I_new[i*ns])
## Escribir en la terminal %matplotlib auto
fig=plt.figure()
ax=plt.gca()
def actualizar(i):
    ax.clear()
    plt.plot(x,I_new[i*10],"blue")
    plt.title("I(x,t) versus x, t=0,6,...,60.")
    plt.xlim(-3,3)
    plt.ylim(0,1.5)
ani=animation.FuncAnimation(fig,actualizar,range(int(len(I_new)/10)))
plt.show()
#%%Se grafica la solucion de la poblacion de Recuperados
for i in range(nplot):
    plt.grid(True)
    plt.title("R(x,t) versus x, t=0,6,...,60.")
    plt.xlabel('x')
    plt.ylabel("R(x,t) t=0,6,...,60.")
    plt.plot(x,R_new[i*ns])
## Escribir en la terminal %matplotlib auto
fig=plt.figure()
ax=plt.gca()
def actualizar(i):
    ax.clear()
    plt.plot(x,R_new[i*10],"blue")
    plt.title("R(x,t) versus x, t=0,6,...,60.")
    plt.xlim(-3,3)
    plt.ylim(0,15)
ani=animation.FuncAnimation(fig,actualizar,range(int(len(R_new)/10)))
plt.show()