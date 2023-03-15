import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from dwave.system import DWaveSampler, EmbeddingComposite
import dwave.inspector as inspector
from dimod import BinaryQuadraticModel as BQM

# BQM representa un modelo cuadratico binario (Ising o QUBO). podemos pasar de qubo a ising como se ve.
# Se puede definir indicando cada vertice y arista con su respectivo peso.
bqm = BQM({0: 7.6, 4: 3.6},{(0, 4): -5}, 0.0, 'BINARY')
print(bqm.to_ising())

# DEFINIMOS UN MANUAL EMBEDDING. PASAMOS LOS qubit_biases Y coupler_strenghts.
# Se mapea en la QPU como: edge_bias = coupler_strenght * 0.25
# Y qubit_bias = sum(edge_bias) + qubit_biases*0.5


# El grafo pasado debe ser un subgrafo del chimera. Lo que es logico ya que al indicar las conexiones
# cuando escribimos por ejemplo (0, 4) estamos refiriendonos a la arista fisica (0, 4) de la QPU.
# Si probamos agregar la arista (4, 5) por ejemplo tira error.
def manual_embedding(qubit_biases, coupler_strengths):

    sampler_manual = DWaveSampler(solver={'topology__type': 'chimera'})
    Q = {**qubit_biases, **coupler_strengths}
    sampleset = sampler_manual.sample_qubo(Q, num_reads=2)
    inspector.show(sampleset)


qubit_biases = {(0, 0): -1, (1, 1): 0, (4, 4): 0, (5, 5): 0}
# One vertex is mapped to the tree {0, 5}
coupler_strengths = {(0, 5): -3, (0, 4): 1, (1, 4): 1, (1, 5): 1}
manual_embedding(qubit_biases, coupler_strengths)



'''
bqm_embedded = BinaryQuadraticModel({1144: 7.6, 1148: 3.6, 1272: 3.6, 1145: 3.4, 1273: 3.4, 1150: 7.6, 1142: 3.6, 1147: 3.6, 1277: 7.6, 1269: 3.6, 1275: 3.6, 1146: 7.6, 1149: 3.6, 1274: 3.6, 1139: 7.6, 1140: 3.6, 1141: 3.6, 1270: 7.6, 1278: 3.6, 1267: 3.6, 1143: 7.6, 1151: 3.6, 1137: 3.6, 1271: 7.6, 1265: 3.6, 1279: 3.6, 1136: 3.4, 1264: 3.4}, {(1144, 1150): 0.2, (1144, 1149): 0.2, (1144, 1272): -8, (1144, 1148): -8, (1148, 1145): 0.4, (1148, 1147): 0.2, (1148, 1146): 0.2, (1148, 1140): 0.4, (1272, 1277): 0.4, (1272, 1278): 0.4, (1145, 1150): 0.4, (1145, 1149): 0.4, (1145, 1151): 0.4, (1145, 1273): -8, (1273, 1277): 0.4, (1273, 1279): 0.4, (1150, 1142): -8, (1150, 1147): -8, (1142, 1139): 0.4, (1142, 1137): 0.2, (1142, 1136): 0.4, (1147, 1275): 0.4, (1147, 1151): 0.2, (1277, 1269): -8, (1277, 1275): -8, (1269, 1267): 0.2, (1269, 1265): 0.2, (1269, 1264): 0.4, (1275, 1278): 0.2, (1275, 1279): 0.2, (1146, 1151): 0.4, (1146, 1274): -8, (1146, 1149): -8, (1149, 1141): 0.4, (1274, 1278): 0.4, (1274, 1279): 0.4, (1139, 1267): 0.4, (1139, 1143): 0.13333333333333333, (1139, 1140): -8, (1139, 1141): -8, (1140, 1137): 0.13333333333333333, (1140, 1136): 0.2, (1141, 1137): 0.13333333333333333, (1141, 1136): 0.2, (1270, 1265): 0.2, (1270, 1264): 0.4, (1270, 1267): -8, (1270, 1278): -8, (1267, 1271): 0.2, (1143, 1136): 0.4, (1143, 1137): -8, (1143, 1151): -8, (1137, 1265): 0.4, (1271, 1264): 0.4, (1271, 1265): -8, (1271, 1279): -8, (1136, 1264): -8}, 0.0, Vartype.BINARY)

sampler = DWaveSampler(solver={'lower_noise': False, 'qpu': True})
raw_result = sampler.sample(bqm_embedded,num_reads=20,annealing_time=1,auto_scale=False)
'''