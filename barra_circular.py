import numpy as np
import sympy as sym

class Matrix:
    def __init__(self, n_var, list_ec):
        self._n_var = n_var
        self._list_ec = list_ec
        
    def T_symb(self):
        T = []
        
        for i in range(self._num_var):
            i = str(i)
            Ti = sym.symbols('T' + i)
            T.append(Ti)
            return T
        
    def eqs_symb(self):
        eqs = []
        
        for j in range(self.num_var):
            eqs.append(self._list_ec[j])
            return eqs

# Propiedades de la barra
material = 'acero'
k = 63.9 # conductividad térmica en W/m.ºC
alpha = 18.8E-6 # difusividad térmica en m^2/s
h = 25 # coeficiente de convección en W/m^2.ºC
q_dot = 5000 # flujo de calor W/m^2
L = 1 # longitud de la barra en m
radio = 0.02 # radio de la barra en m
delta_r = 0.001 # distancia entre nodos en m
delta_z = delta_r # nodos equiespaciados
M = int(L // delta_z) # cantidad de secciones en z
N = int(radio // delta_r) # cantidad de secciones en r

# Parte temporal
t = 600 # tiempo total en s
delta_t = 10 # intervalo de tiempo
time_iter = t // delta_t # iteraciones temporales

T_amb = 20 # en ºC

tau = (delta_r)**2 / (alpha * delta_t)


for k in range(0, time_iter, 1):
    for m in range(0, radio, delta_r):
        for n in range(0, L, delta_z):
            eq1 = -tau * t[k+1][m-1][n] - tau * t[k+1][m+1][n] - tau * t[k+1][m][n-1] - tau * t[k+1][m][n+1] + (4 * tau +1) * t[k+1][m][n] - t[k][m][n] 


     
