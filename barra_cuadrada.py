import numpy as np
import matplotlib.pyplot as plt

# Propiedades de la barra
material = 'acero'
k = 63.9 # conductividad térmica en W/m.ºC
alpha = 18.8E-6 # difusividad térmica en m^2/s
h = 25 # coeficiente de convección en W/m^2.ºC
L = 1 # longitud de la barra en m
width = 0.05 # ancho de la barra en m
delta_x = 0.025 # distancia entre nodos en m
delta_y = delta_x # nodos equiespaciados
M = int(L // delta_x) # cantidad de secciones en x
N = int(L // delta_y) # cantidad de secciones en y

# Parte temporal
t = 600 # tiempo total en s
delta_t = 10 # intervalo de tiempo
time_iter = t // delta_t # iteraciones temporales

# Condiciones iniciales
T_initial = 20.0 # temperatura inicial en ºC

# Condiciones de contorno en ºC
T_top = 20.0
T_bottom = 20.0
T_left = 150.0
T_right = 20.0

# Red nodal
dot_x = M + 1
dot_y = N + 1

T_t = np.empty((dot_x, dot_y)) # armamos un arreglo igual a la red nodal
T_t.fill(T_initial) # ponemos todos los nodos a la temp inicial

# Ponemos las condiciones de contorno en el arreglo
T_t[:, (L-1)] = T_left



#Definimos la clase
class Matrix:
    _M = 0
    _N = 0
    _elem = None

#vamos a inicializar los atributos del objeto
    def __init__(self, M, N): 
        self._M = M
        self._N = N
        self._elem = []
        for i in range(self._M):
            self._elem.append([])
            for j in range(self.n):
                self._elem[i].append(0)
                
        print(self._elem)
        