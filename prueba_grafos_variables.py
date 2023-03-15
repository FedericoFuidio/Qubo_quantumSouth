# Probamos meterle a D-wave un problema con la ecuacion a resolver:

from dwave.system.samplers import DWaveSampler

from dwave.system.composites import EmbeddingComposite

sampler = EmbeddingComposite(DWaveSampler(endpoint='', token='', solver=''))

query3 = {('x1', 'x2'): 2, ('x1', 'z'): -2, ('x2', 'z'): -2, ('z', 'z'): 2, ('z', 'a'): 2, ('a', 'a'): -1, ('x1', 'x1'):1, ('x2', 'x2'):1, ('x2', 'd'):-2, ('x1', 'd'):-2, ('d', 'd'):1}

response = sampler.sample_qubo(query3, num_reads=200)

for datum in response.data(['sample', 'energy', 'num_occurrences']):

    print(datum.sample, "Energy: ", datum.energy, "Occurrences: ", datum.num_occurrences)