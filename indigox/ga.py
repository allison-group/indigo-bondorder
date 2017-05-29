from collections import deque
from itertools import combinations
from multiprocessing import Pool
from random import choice, random, sample
from time import perf_counter

from indigox.config import (ELECTRON_PAIRS, MUTATE_PROB, MIN_GENERATIONS,
                            MAX_GENERATIONS, BRUTEFORCE_CUTOFF, CONVERGENCE,
                            POP_SIZE, ELITEISM_SIZE, BREEDING_ELITEISM,
                            SEED_COUNT, TIMEOUT, NUM_PROCESSES)
from indigox.misc import (formal_charge, node_energy, HashBitArray,
                          BondOrderAssignment, graph_to_dist_graph,
                          electron_spots, electrons_to_add, locs_sort,
                          bitarray_to_assignment)
from scipy.misc import comb

import numpy as np


def _ga_calc_energy(G, a, locs):
    for n in G:
        if len(n) == 1:
            G.node[n]['e-'] = 0
        elif len(n) == 2:
            G.node[n]['e-'] = 2
    for i in range(a.length()):
        if a[i] and ELECTRON_PAIRS:
            G.node[locs[i]]['e-'] += 2
        elif a[i]:
            G.node[locs[i]]['e-'] += 1
    # calculate formal charges 
    for n in G:
        if len(n) == 1:
            G.node[n]['fc'] = formal_charge(G, n)
    ene = sum(node_energy(G, n) for n in G)
    return a, round(ene,5)

def _ga_brute_force(args):
    G, combo, locs = args
    a = HashBitArray(len(locs))
    a.setall(False)
    for i in combo:
        a[i] = True
    return _ga_calc_energy(G, a, locs)

def _ga_algroithm(args):
    G, locs, a, b, uniques = args
    c, d = _ga_crossover(a, b, uniques)
    a_ene = _ga_calc_energy(G, a, locs)
    b_ene = _ga_calc_energy(G, b, locs)
    c_ene = _ga_calc_energy(G, c, locs)
    d_ene = _ga_calc_energy(G, d, locs)
    return a_ene, b_ene, c_ene, d_ene
    
def _ga_mutation(bits, uniques):
    
    truths = [i for i in range(bits.length()) if bits[i]]
    drop_bit = choice(truths)
    try:
        pick_bit = choice(list(set(range(bits.length())).difference(truths, uniques[drop_bit])))
    except IndexError:
        # If can't mutate to something unique, just mutate in place
        pick_bit = choice(list(set(range(bits.length())).difference(truths)))
    a = bits.copy()
    a[drop_bit] = False
    a[pick_bit] = True
    return a

def _ga_crossover(a, b, uniques):
    if random() < MUTATE_PROB:
        a = _ga_mutation(a, uniques)
    if random() < MUTATE_PROB:
        b = _ga_mutation(b, uniques)
        
    cross_points = []
    a_count, b_count = 0, 0
    for i in range(a.length()):
        a_count += a[i]
        b_count += b[i]
        if a_count == b_count:
            cross_points.append(i+1)
    try:
        cross = choice(cross_points)
    except IndexError:
        cross = a.length()
    
    c = a[:cross] + b[cross:]
    d = b[:cross] + a[cross:]
        
    return c, d 
    
class GeneticAlogrithm(BondOrderAssignment):
    def __init__(self, G):
        self.init_G = G
                  
    def initialise(self):
        self.tries = 0
        self.success = 0
        self.population = {}
        self.already_seen = {}
        self.gen_count = 0
        self.G = graph_to_dist_graph(self.init_G)
        self.locs = locs_sort(electron_spots(self.init_G), self.G)
        self.target = electrons_to_add(self.init_G)
        self.last_n_energies = deque(maxlen=CONVERGENCE)
        
        self.unique_slots = {}
        for i in range(len(self.locs)):
            self.unique_slots[i] = set()
            for j in range(len(self.locs)):
                if self.locs[i] == self.locs[j]:
                    self.unique_slots[i].add(j)  
                elif j > i:
                    break

    def run(self):
        self.start_time = perf_counter()
        self.initialise()
        self.pool = Pool(processes=NUM_PROCESSES)
        if comb(len(self.locs), self.target) < BRUTEFORCE_CUTOFF:
            self._brute_force()
        else:
            self._initialise()
            while self.gen_count < MAX_GENERATIONS:
                self.gen_count += 1
                new_pop = self.population.copy()
                results = self.pool.imap_unordered(_ga_algroithm,
                                    ((self.G, self.locs, a, b, self.unique_slots) for a, b in self._selections()),
                                    chunksize=25)
                for group in results:
                    if perf_counter() - self.start_time > TIMEOUT:
                        self.pool.terminate()
                        raise TimeoutError
                    for a, ene in group:
                        new_pop[a] = ene
                        
                self.population = self._killoff(new_pop)
                for k, v in self.population.items():
                    self.already_seen[k] = v
                current_min = int(round(min(self.population.values()), 5)) * 1e5
                self.last_n_energies.append(current_min)
                if (self.gen_count >= MIN_GENERATIONS 
                    and len(self.last_n_energies) == CONVERGENCE
                    and len(set(self.last_n_energies)) == 1):
                    break
        self.pool.terminate()
        opt = min(self.already_seen, key=lambda x:self.already_seen[x])
        bitarray_to_assignment(self.init_G, opt, self.locs)
        return self.init_G, self.already_seen[opt]
      
    def _killoff(self, pop):
        alive = {}
        while len(alive) < POP_SIZE:
            if perf_counter() - self.start_time > TIMEOUT:
                self.pool.terminate()
                raise TimeoutError
            sorted_pop = sorted(x for x in pop if x not in alive)
            z = sum(np.exp(-pop[x]) for x in sorted_pop)
            select_probs = [self._xi_prob(pop[x], z, ELITEISM_SIZE) for x in sorted_pop]
            for i in range(1,len(select_probs)):
                select_probs[i] = select_probs[i] + select_probs[i-1]
            select_probs[-1] = 1
            r1 = random()
            i = next(idx for idx, v in enumerate(select_probs) if v >= r1)
            k = sorted_pop.pop(i)
            alive[k] = pop[k]
        return alive 
      
    def _xi_prob(self, ene, z, XI):
        return (np.exp(-ene) / z) * XI + (1 - XI) / POP_SIZE
    
    def _selections(self):
        sorted_pop = sorted(self.population)
        z = sum(np.exp(-self.population[x]) for x in sorted_pop)
        select_probs = [self._xi_prob(self.population[x], z, BREEDING_ELITEISM) for x in sorted_pop]
        for i in range(1,len(select_probs)):
            select_probs[i] = select_probs[i] + select_probs[i-1]
        select_probs[-1] = 1
        seletions = []
        for _ in range(POP_SIZE):
            r1 = random()
            i = sorted_pop[next(idx for idx, v in enumerate(select_probs) if v >= r1)]
            r2 = random()
            j = sorted_pop[next(idx for idx, v in enumerate(select_probs) if v >= r2)]
            while j == i:
                if perf_counter() - self.start_time > TIMEOUT:
                    self.pool.terminate()
                    raise TimeoutError
                r2 = random()
                j = sorted_pop[next(idx for idx, v in enumerate(select_probs) if v >= r2)]
            seletions.append((i,j))
        return seletions
    
    def _brute_force(self):
        self.log.warning('Running GeneticAlgoritm with brute force')
        results = self.pool.imap_unordered(_ga_brute_force,
                                ((self.G, x, self.locs) for x in combinations(range(len(self.locs)), self.target)))
        for a, ene in results:
            if perf_counter() - self.start_time > TIMEOUT:
                self.pool.terminate()
                raise TimeoutError
            self.already_seen[a] = ene
        
    def _initialise(self):
        if SEED_COUNT:
            self._seeding()
        while len(self.population) < POP_SIZE:
            a = HashBitArray(len(self.locs))
            if perf_counter() - self.start_time > TIMEOUT:
                self.pool.terminate()
                raise TimeoutError
            while True:
                if perf_counter() - self.start_time > TIMEOUT:
                    self.pool.terminate()
                    raise TimeoutError
                a.setall(False)
                locs = sample(range(len(self.locs)), k=self.target)
                for i in locs:
                    a[i] = True
                if a not in self.population:
                    ene = _ga_calc_energy(self.G, a, self.locs)[1]
                    self.population[a] = ene
                    self.already_seen[a] = ene
                    break
    
    def _seeding(self):
        seed = HashBitArray(len(self.locs))
        seed.setall(False)
        base_energy = _ga_calc_energy(self.G, seed, self.locs)[1]
        all_locs = list(range(len(self.locs)))
        while seed.count() < self.target:
            if perf_counter() - self.start_time > TIMEOUT:
                self.pool.terminate()
                raise TimeoutError
            energy_diffs = {}
            for i in all_locs:
                seed[i] = True
                energy_diffs[i] = (_ga_calc_energy(self.G, seed, self.locs)[1]
                                   - base_energy)
                seed[i] = False
            min_i = min(energy_diffs, key=lambda x: energy_diffs[x])
            seed[min_i] = True
            base_energy += energy_diffs[min_i]
            all_locs.remove(min_i)
            
        self.population[seed] = base_energy
        self.already_seen[seed] = base_energy
        max_muts = 5
        while len(self.population) < SEED_COUNT:
            seed_x = _ga_mutation(seed, self.unique_slots)
            if perf_counter() - self.start_time > TIMEOUT:
                self.pool.terminate()
                raise TimeoutError
            if seed_x not in self.population:
                ene_x = _ga_calc_energy(self.G, seed_x, self.locs)[1]
                self.population[seed_x] = ene_x
                self.already_seen[seed_x] = ene_x
                max_muts -= 1
                if ene_x < base_energy or not max_muts:
                    seed = seed_x
                    base_energy = ene_x
                    max_muts = 5
        
            
