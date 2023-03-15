import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector as inspector
from dimod import BinaryQuadraticModel as BQM

# p = 4, n = 16

def chimera_column(k):
    return int(np.floor(k/8) % 16)

def chimera_row(k):
    return int(np.floor(k/(8*16)))

def qubit_column(k):
    return int(4*chimera_column(k) + (k % 4))

def qubit_row(k):
    return int(4*chimera_row(k) + (k % 4))


# La variable n del problema se mapea a los qbits diagonales
def qubits_variable(n):
    chimera = int(np.floor(n/4)) # Nos dice el numero de la chimera diagnoal (chimera, chimera)
    qubit_chimera = int(n) % 4 # Vamos a usar qubit_chimera y qubit_chimera + 4 como diagonal
    qubit_number = chimera * (8*16) + 8*(chimera % 8)
    return qubit_number + qubit_chimera


def add_horizontal_edges(qubit, edges, n_mod4, chains):
    # edges: del 1 al 4

    for i in range(0, edges):
        chains[(qubit - n_mod4 - (4 - i), qubit)] = 1
    
    return chains

def add_vertical_edges(qubit, edges, n_mod4, chains):
    # edges: del 1 al 4

    for i in range(0, edges):
        chains[(qubit, qubit - n_mod4 + (4 + i))] = 1
    
    return chains
    

# Retorna todos los qubits relacionados con la variable n.
def qubits_dict(n, variables): 
    # n ---> numero de variable en el problema (x_0, x_1, ..., x_n etc)
    # variables ---> Maximo valor de los indices de variables. arrancando de cero

    chimera = int(np.floor(n/4)) # Chimera diagona donde se unen la cadena horizonal y vertical
    dimentions = int(np.floor(variables/4)) # Chimeras diagonales a usar (indice de maxima chimera diagonal)
    max_qubit = qubits_variable(variables) # Qubit de mayor valor a usar:


    # n debe ser a lo sumo igual a la cantidad de variables totales
    if(n > variables):
        return

    n_mod4 = n % 4 # Vamos a usar varias veces la variable, esta bueno tenerla.
    # Qubits de inicio
    Vhori = chimera*(8*16) + n_mod4 + 4
    Vvert = chimera*(8) + n_mod4

    print(Vvert, Vhori)

    qubits_hori = []
    qubits_vert = []

    qubits = dict()
    chains = dict()
    for i in range(0, dimentions + 1):

        qubits_hori.append(Vhori + (8*i)) 
        qubits_vert.append(Vvert + (8*16)*i)
        qubits[(qubits_hori[i], qubits_hori[i])] = 1
        qubits[(qubits_vert[i], qubits_vert[i])] = 1

        # Debemos agregarle los pesos
        
        if (i > chimera):
            chains = add_horizontal_edges(qubits_hori[i], 4, n_mod4, chains)
            chains = add_vertical_edges(qubits_vert[i], 4, n_mod4, chains)

        if (i > 0):
            chains[(qubits_hori[i], qubits_hori[i-1])] = 1
            chains[(qubits_vert[i], qubits_vert[i-1])] = 1



    # conectamos tambien los elementos Floor(n/4) (los que estan en la misma chimera)

    chains[(qubits_vert[chimera], qubits_hori[chimera])] = 1

    return [qubits, chains]
