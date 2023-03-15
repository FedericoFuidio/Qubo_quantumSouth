import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector as inspector

# PONEMOS UNA MATRIZ Y VEMOS COMO CORRE MANUALMENTE
# ESPECIFICAMOS chain_strenght en embedding.sample_qubo

identity = np.identity(5)
identity[0][0] = -1

matriz = np.matrix([[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [1, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
                    [1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0],
                    [0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1],
                    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]])

matriz = np.matrix([[0, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1],
                    [8, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1],
                    [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                    [1, 0, 0, 0, 1, -10, 1, 1, 1, 1, 1],
                    [1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1],
                    [1, 0, 0, -10, 1, 0, 1, 1, 1, 1, 1],
                    [1, 0, 0, 1, 1, 1, 0, 1, 1, 1, 1],
                    [1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1],
                    [1, 0, 0, 1, 1, 1, 1, 1, 0, 1, 1],
                    [1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 1],
                    [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]])

# Probamos con un grafo completo de (16*4)*(16*4) (maximo posible. Si hay vertices rotos no lo hara)
matriz = np.ones((16 + 1, 16 + 1))

sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})
embedding = EmbeddingComposite(sampler_manual)
running = embedding.sample_qubo(matriz, num_reads=2, chain_strength=0.5)
inspector.show(running)


'''
sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})

# Check if quibits and couplers are available: 
qubits_free = all(qubit in sampler_manual.nodelist for qubit in [0, 1, 4, 5])
couplers_free = all(coupler in sampler_manual.edgelist for coupler in [(0, 4), (0, 5), (1, 4), (1, 5)])



qubit_biases = {(0, 0): 0, (1, 1): 0, (2, 2): 0, (3, 3): 0, (4, 4): 0, (5, 5): 0, (6, 6): 0, (7, 7): 0,
             (8, 8): 0, (9, 9): 0, (12, 12): 0}
coupler_strengths = {(0, 4): -1, (0, 6): 1, (1, 5): 1, (1, 6): 1, (2, 4): 1, (2, 6): 1, (3, 6): 1, 
                    (3, 7): 1, (4, 12):1,(8, 12): 1, (9, 12): 1}
                    
Q = {**qubit_biases, **coupler_strengths}

#sampleset = sampler_manual.sample_qubo(Q, num_reads=2)
#inspector.show(sampleset)
'''