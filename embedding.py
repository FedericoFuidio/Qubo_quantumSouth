import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector as inspector
from mapping_qbits import *

def Merge(dict_1, dict_2):
	result = dict_1 | dict_2
	return result

# Embedding a K_{3} into the chimera manually and authomaticaly:

def auto_embedding(matrix): 

    sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})
    embedding = EmbeddingComposite(sampler_manual)
    running = embedding.sample_qubo(matrix, num_reads=20)
    print(running)
    inspector.show(running)

def manual_embedding(qubit_biases, coupler_strengths):

    sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})
    Q = {**qubit_biases, **coupler_strengths}
    sampleset = sampler_manual.sample_qubo(Q, num_reads=2)
    print(sampleset)
    inspector.show(sampleset)

map = qubits_dict(0, 4)
qubits = map[0]
chains = map[1]

map2 = qubits_dict(1, 4)
qubits2 = map2[0]
chains2 = map2[1]

map3 = qubits_dict(2, 4)
qubits3 = map3[0]
chains3 = map3[1]

map4 = qubits_dict(3, 4)
qubits4 = map4[0]
chains4 = map4[1]

map5 = qubits_dict(4, 4)
qubits5 = map5[0]
chains5 = map5[1]

Q = Merge(qubits, qubits2)
C = Merge(chains, chains2)
Q = Merge(Q, qubits3)
C = Merge(C, chains3)
Q = Merge(Q, qubits4)
C = Merge(C, chains4)
Q = Merge(Q, qubits5)
C = Merge(C, chains5)
manual_embedding(Q, C)


