import numpy as np
import pandas as pd

def get_Q_objetivo(profits, num_slack):
  diagonal = np.r_[profits, np.zeros(num_slack)] # Creamos la diagonal de la matriz
  qubo = np.diag(diagonal) # Creamos la matriz diagonal
  return qubo

def get_Q_Peso(max_weight, weights, num_slack):
 
    pesosSlacks = (2 ** (np.arange(num_slack - 1))).astype(float)
    pesosSlacks = np.r_[pesosSlacks, max_weight - sum(pesosSlacks)] # revisar esto.
    #pesosSlacks = np.r_[[1,2,4,8,16],[2]]
    vect = np.r_[weights, pesosSlacks] #Agregamos un menos weights, a ver que sale
    # Por que el -2 * max_weight * np.diag(vect) en la cuenta?
    qubo = np.outer(vect, vect) - 2 * max_weight * np.diag(vect)
    return qubo

def addOne(indice, set):

    set[indice][1] += 1
    return

# Retorna cada elemento en sampleSet, junto con la correspondiente
# cantidad de veces que aparece:
def histograma(sampleSet):

    set = []

    samp = np.array([element[0] for element in sampleSet])
    samp = [samp[i, :] for i in range(0, len(sampleSet))]
    
    for sample in samp:
        print(sample)
        if sample in np.array(row[0] for row in set):
          indice = np.where(np.array(row[0] for row in set) == sample)[0][0]
          addOne(indice, set)

        else:
          set.append([sample, 1])

    return set


