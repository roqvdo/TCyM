import numpy as np
import sympy as sym
from sympy import *
from sympy.matrices import Matrix 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

#class Matrix:
    #def __init__(self, n_var, list_ec):
       # self._n_var = n_var
       # self._list_ec = list_ec
        
    #def T_symb(self):
       #T = []
        
       #for i in range(self._num_var):
         # i = str(i)
         # Ti = sym.symbols('T' + i)
          #T.append(Ti)
          #return T
        
    #def eqs_symb(self):
     #   eqs = []
        
       # for j in range(self.num_var):
       #     eqs.append(self._list_ec[j])
        #    return eqs

# Propiedades de la barra
material = 'acero'
k = 63.9 # conductividad térmica en W/m.ºC
alpha = 18.8E-6 # difusividad térmica en m^2/s
h = 25 # coeficiente de convección en W/m^2.ºC
q_dot = 5000 # flujo de calor W/m^2
L = 1 # longitud de la barra en m
radio = 0.02 # radio de la barra en m
delta_r = 0.005 # distancia entre nodos en m
delta_z = delta_r # nodos equiespaciados
M = int(L // delta_z) # cantidad de secciones en z
N = int(radio // delta_r) # cantidad de secciones en r

# Parte temporal
t = 600 # tiempo total en s
delta_t = 1 # intervalo de tiempo
time_iter = t // delta_t # iteraciones temporales

T_amb = 20 # en ºC

tau = (alpha * delta_t) / (delta_r)**2 

# Red nodal para la temperatura inicial
t = np.zeros((M + 1, N + 1))
for m in range(M+1):
  for n in range(N+1):
    t[m][n] = T_amb
    

# Definición de las variables simbólicas
T = []
varibles = (M + 1) * (N + 1)

for i in range(varibles):
  i = str(i)
  Ti = sym.symbols('T' + i)
  T.append(Ti)
  
# Red nodal de temperaturas
g= Matrix(T).reshape(M+1,N+1)

# Ecuaciones para los puntos de la red nodal
ecuaciones = []

for m in range(M + 1):
    for n in range(N + 1):
        if 0 < m < M and 0 < n <N: # Nodos interiores
            E1 = - tau * g[m-1, n] - tau * g[m+1, n]- tau * g[m, n-1]- tau * g[m, n+1] + (4 * tau + 1) * g[m, n] - t[m][n]
            ecuaciones.append(E1)
        elif 0 < m < M and n == N: # Nodos frontera superiores
            E2 = ((2 * tau * h * delta_r) / k) + (4 * tau + 1) * g[m, n] - 2 * tau * g[m, n-1] - tau * g[m-1, n] - tau * g[m+1, n] - ((2 * tau * h * delta_r) / k) * T_amb - t[m][n]
            ecuaciones.append(E2)
        elif 0 < m < M and n == 0: # Nodos frontera inferiores
            E3 = - tau * g[m-1, n] + (4 * tau + 1) * g[m, n] - tau * g[m+1, n] - 2 * tau * g[m, n+1] - t[m][n]
            ecuaciones.append(E3)
        elif m == 0 and 0 < n < N: # Nodos frontera izquierda
            E4 = - 2 * tau * g[m+1, n] + (4 * tau + 1) * g[m, n]- tau * g[m, n+1] - ((2 * tau * q_dot * delta_r) / k) - t[m][n]
            ecuaciones.append(E4)
        elif m == 0 and n == 0: # Nodo inferior izquierdo
            E5 = - 2 * tau * g[m+1, n] + (4 * tau + 1) * g[m,n] - 2 * tau * g[m,n+1] - ((2 * tau * q_dot * delta_r) / k) - t[m][n]
            ecuaciones.append(E5)
        elif m == M and 0 < n < N: # Nodos frontera derecha
            E6 = (((2 * tau * h * delta_r) / k) + 4 * tau + 1) * g[m, n] - 2 * tau * g[m-1, n] - tau * g[m, n-1] - ((2 * tau * h * delta_r) / k) * T_amb - t[m][n]
            ecuaciones.append(E6)
        elif m == M and n == 0: # Nodo inferior derecho
            E7 = (((2 * tau * h * delta_r) / k) + 4 * tau + 1) * g[m, n] - 2 * tau * g[m-1, n] - tau * g[m, n+1] - ((2 * tau * h * delta_r) / k) * T_amb - t[m][n]
            ecuaciones.append(E7)
        elif m == 0 and n == N: # Nodo superior izquierdo
            E8 = - 4 * tau * g[m+1,n] - 4 * tau * g[m, n-1] + (((2 * tau * h * delta_r) / k) + 8 * tau + 1) * g[m, n] - ((2 * tau * h * delta_r) / k) * T_amb - ((2 * tau * q_dot * delta_r) / k) - t[m][n]
            ecuaciones.append(E8)
        elif m == M and n == N: # Nodo superior derecho
            E9 = (((4 * tau * h * delta_r) / k) + 8 * tau + 1) * g[m, n] - 4 * tau * g[m-1, n] - 4 * tau * g[m, n-1] - ((4 * tau * h * delta_r) / k) * T_amb - t[m][n]
            ecuaciones.append(E9)
            
# Definimos la iteración para las temperaturas        
def temp_iteration(time_iter):
  for k in range(time_iter):
    for m in range(M+1):
      for n in range(N+1):
        if k == 0:
          t[m][n] = T_amb
          A, b = sym.linear_eq_to_matrix(ecuaciones, T)
          system = (A, b)
          T_k, = sym.linsolve(system, T) # la coma es para poder acceder a los elementos
          solution = Matrix(T_k).reshape(M+1,N+1)
        else:
          t[m][n] = solution[m][n]
          A, b = sym.linear_eq_to_matrix(ecuaciones, T)
          system = (A, b)
          T_k, = sym.linsolve(system, T)
          solution = Matrix(T_k).reshape(M+1,N+1)
          
  return solution
        
    
def plotheatmap(solution, time_iter):
  plt.clf()
  plt.title(f"Temperatura para t = {k*delta_t:.3f} segundos")
  plt.xlabel("z")
  plt.ylabel("r")
  
  # This is to plot solution (solution at time-step k)
  plt.pcolormesh(solution, cmap=plt.cm.jet, vmin=0, vmax=100)
  plt.colorbar()
  
  return plt

solution = temp_iteration(time_iter)

def animate(k):
    plotheatmap(solution[k], k)

anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=time_iter, repeat=False)
anim.save("solucion_barra_cilindrica.gif")

     
