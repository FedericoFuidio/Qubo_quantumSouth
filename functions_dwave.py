import os
import time
import numpy as np
import pandas as pd
import math
import dimod
import warnings
warnings.filterwarnings('ignore')

def get_Q_objetivo(profits, num_slack):
  diagonal = np.r_[profits, np.zeros(num_slack)] # Creamos la diagonal de la matriz
  qubo = np.diag(diagonal) # Creamos la matriz diagonal
  return qubo

def get_Q_Peso(max_weight, weights, num_slack_weights, num_slack):
 
    pesosSlacks = (2 ** (np.arange(num_slack_weights - 1))).astype(float)
    pesosSlacks = np.r_[pesosSlacks, max_weight - sum(pesosSlacks)] # revisar esto.
    # pesosSlacks = np.r_[[1,2,4,8,16],[2]]
    vect = np.r_[weights, pesosSlacks, np.zeros(num_slack - num_slack_weights)] #Agregamos un menos weights, a ver que sale
    
    qubo = np.outer(vect, vect) - 2 * max_weight * np.diag(vect)
    return qubo

def get_Q_Volumen(max_volume, volumes, num_slack_volumes, num_slack):
    volumenSlacks = (2 ** (np.arange(num_slack_volumes - 1))).astype(float)
    volumenSlacks = np.r_[volumenSlacks, max_volume - sum(volumenSlacks)]
    vect = np.r_[volumes, np.zeros(num_slack - num_slack_volumes), volumenSlacks]

    qubo = np.outer(vect, vect) - 2 * max_volume * np.diag(vect)
    return qubo

# Vector menor energia (Tendremos que encontrarlo entre todas las soluciones de Dwave):
def lowest_energy(sampleset):  # Finds the lowest energy solution
    
    # Description: given a full sampleset (tuples of the form (solution, energy) finds the lowest energy SAMPLE.
    # INPUTS:
    # Sampleset: a sampleset of the form list((solution, energy))

    # OUTPUTS:
    # best: a tuple of the form (solution, energy)

    if len(sampleset):
        #energies = np.array(sampleset)[:, 1]  # energias
        energies = [row[1] for row in sampleset]
        index = np.argmin(energies)  # indice de la de menor energia
        ret = sampleset[index]  # solucion de menor energia
        return ret

    else:
        return None


# Funciones de Checkeo:
def wrong_slack(sample, max_weight, weights, num_slack_weights, num_slack):
    
    # INPUTS:
    # Sample (list): vector x que queremos checkear
    # max_weight (int): peso máximo
    # weights (list): lista de pesos
    # num_slack (int): cantidad de slacks para representar max_weight
    
    # OUTPUTS:
    # False si la solucion es incorrecta
 
    ret = False  # valid solution if true

    # creamos el vector "vect" de slacks teorico, mencionado en el cuaderno de jupyter "qubo.ipynb"

    pesosSlacks = 2 ** (np.arange(num_slack_weights - 1))
    pesosSlacks = np.r_[pesosSlacks, max_weight - sum(pesosSlacks)]
    vect = np.r_[-weights, pesosSlacks, np.zeros(num_slack - num_slack_weights)] # nos sirve para checkear ACA TAMBIEN CAMBIO SIGNO
    solWeightSlack = sample * vect 


    ws = sum(solWeightSlack) # en ws queda guardado la suma de los pesos cargados + el numero que representan las slacks
    print("ws = ", ws)
    #if abs(ws - max_weight) > 0.001
    if abs(ws) > 0.001: # si la diferencia es casi 0... #DICE MAYOR EN CODIGO ORIGINAL
        ret = True # las slacks se prendieron correctamente
    return ret


def check_weight(sample, max_weight, weights, num_slack_weights, num_slack):
    # Description:
    # 

    # INPUTS:
    # Sample (list): vector x que queremos checkear
    # max_weight (int): peso máximo
    # weights (list): lista de pesos
    # num_slack (int): cantidad de slacks para representar max_weight
    # OUPUTS:
    # True if solution checks weights inequality, False otherwise

    ret = True
    error_slack = wrong_slack(sample, max_weight, weights, num_slack_weights, num_slack)
    sample = sample[0:len(sample) - (num_slack)]  # Me quedo con la solucion, ya no me sirven las slacks

    loaded_weight = sum(sample * weights)

    if error_slack:  # si hay error en las slacks, ya descarto de una
        print("weight slack  error")
        ret = False
    else:  # si no hay error en las slacks, me fijo el peso cargado
        if loaded_weight > loaded_weight:
            ret = False
            print("weigh exceeded")

    return ret


#Funcion de simulacion (envia a Dwave): (¿Hacer un sort?)
def sendToDwave(qubo, shots=100, aggregate = True):
    # Description: functions that solves a particular qubo problem

    # INPUT:
    # qubo: (matrix) representation of the xt*Q*x problem
    # shots: (int) cantidad de simulaciones a realizar
    # aggregate: (boolean) indica si queremos que solo retorne soluciones distintas
    # sort: (boolean) indica si queremos al ordenado por mayor cantidad de apariciones

    # OUTPUT:
    # sampleset: array of tuples of the form (solution, energy, num_occurrences) of length "shots" containing posible (but not neccesarily feasible) solutions
    # or (solution, energy) if aggregate is False (aggregate is True by default)
    
    tic = time.perf_counter() # for time measuring
    sampleset = dimod.SimulatedAnnealingSampler().sample_qubo(qubo, num_reads=shots)
    if(aggregate): # Solo tenemos en cuenta soluciones DISTINTAS, retornamos tambien la cantidad de ocurrencias de cada solucion
        sampleset = sampleset.aggregate() # solo agrega soluciones DIFERENTES. 
        sampleset = [(sample, energy, cant) for sample, energy, cant in zip(sampleset.record.sample, sampleset.record.energy, sampleset.record.num_occurrences)]
    
    else: # Retornamos todo, no pasamos el numero de ocurrencias
        sampleset = [(sample, energy) for sample, energy in zip(sampleset.record.sample, sampleset.record.energy)]

    toc = time.perf_counter() # for time measuring
    print(f"Simmulating {shots} instances of annealing took: {(toc-tic)}s")
    # print("Sampleset sin filtrar: ", sampleset)
    return sampleset


#Funcion de filtrado:
def filterKnapSampleset(sampleset,  max_weight, weights,  max_volume, volumes, num_slack_weights, num_slack):
    num_slack_volumes = num_slack - num_slack_weights
    # Description:
    # given a FULL SAMPLESET for knapsack problem, filters the invalid samplesets

    # INPUTS:
    # sampleset: raw sampleset (list of tuples) returned from dwave sampler, the structure is: (solution, energy)
    # max_weight: (int) maximum weight
    # weights (list): lista de pesos
    # num_slack (int): cantidad de slacks para representar max_weight
    

    # OUTPUTS:
    # feasibleSamples: sampleset (list tuples) with the valid solutions, the structure is: (validSolution, energy)

    feasibleSamples = []
    cantidadValidas = 0

    i = 0
    for sample, energy in sampleset:
        print("-----------------------")
        print("Checking solution: ", i)

        weightRespected = check_weight(sample, max_weight, weights, num_slack_weights, num_slack)  # weightFlag

        if weightRespected:  # if sample is valid:
            print("valid solution")
            feasibleSamples.append((sample, energy))
            cantidadValidas = cantidadValidas + 1
        i = i + 1

    return feasibleSamples



