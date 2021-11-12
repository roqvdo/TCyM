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

# Red nodal de temperaturas
for k in range(time_iter):
    for m in range(M + 1):
        for n in range(N + 1):
            T = np.zeros((M + 1, N + 1))
            

# Ecuaciones
    # Nodos interiores
eq1 = - tau * T[m-1][n][k+1] - tau * T[m+1][n][k+1]- tau * T[m][n-1][k+1]- tau * T[m][n+1][k+1] + (4 * tau + 1) * T[m][n][k+1] - T[m][n][k]
    # Nodos frontera superiores
eq2 = (((2 * tau * h * delta_r) / k) + (4 * tau + 1) * T[m][n][k+1] - 2 * tau * T[m][n-1][k+1] - tau * T[m-1][n][k+1] - tau * T[m+1][n][k+1] - ((2 * tau * h * delta_r) / k) * T_amb - T[m][n][k]
    # Nodos frontera inferiores
eq3 = - tau * T[m-1][n][k+1] + (4 * tau + 1) * T[m][n][k+1] - tau * T[m+1][n][k+1] - 2 * tau * T[m][n+1][k+1] - T[m][n][k]
    # Nodos frontera izquiera
eq4 = - 2 * tau * T[m+1][n][k] + (4 * tau + 1) * T[m][n][k+1]- tau * T[m][n+1][k+1] - ((2 * tau * q_dot * delta_r) / k) - T[m][n][k]
    # Nodo inferior izquierdo
eq5 = - 2 * tau * T[m+1][n][k] + (4 * tau + 1) * T[m][n][k+1] - 2 * tau * T[m][n+1][k+1] - ((2 * tau * q_dot * delta_r) / k) - T[m][n][k]
    # Nodo



# Creamos las variables simbólicas
#num_var = (M + 1) * (N + 1)

#T = []

#for i in range(1, num_var):
    #i = str(i)
    #Ti = sym.symbols('T' + i)
    #T.append(Ti)






# Red nodal de temperaturas
#for k in range(0, time_iter):
   # grid = np.zeros((M + 1, N + 1))
    
    #if k == 0:
     #   for m in range(0, M + 1):
      #      for n in range(0, N + 1):
       #         grid[m][n] = T_amb
    #else:
        
        
    
        


     
