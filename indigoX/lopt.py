from multiprocessing import Pool

from indigoX.config import INIT_WITH_GA, NUM_PROCESSES
from indigoX.exception import IndigoUnfeasibleComputation

from indigoX.misc import (BondOrderAssignment, graph_to_dist_graph, electron_spots,
                   electrons_to_add, locs_sort, HashBitArray, graph_setup,
                   node_energy, bitarray_to_assignment)


class LocalOptimisation(BondOrderAssignment):
    def __init__(self, G):
        self.init_G = G
        
    def initialise(self):
        self.G = graph_to_dist_graph(self.init_G)
        self.target = electrons_to_add(self.init_G)
        self.locs = locs_sort(electron_spots(self.init_G), self.G)
        if INIT_WITH_GA:
            self.init_a = HashBitArray(len(self.locs))
            self.init_a.setall(False)
            base_energy = self.calc_energy(self.init_a)[1]
            all_locs = list(range(len(self.locs)))
            while self.init_a.count() < self.target:
                energy_diffs = {}
                for i in all_locs:
                    self.init_a[i] = True
                    energy_diffs[i] = (self.calc_energy(self.init_a)[1]
                                       - base_energy)
                    self.init_a[i] = False
                min_i = min(energy_diffs, key=lambda x: energy_diffs[x])
                self.init_a[min_i] = True
                base_energy += energy_diffs[min_i]
                all_locs.remove(min_i)
        else:
            self.init_a = HashBitArray(len(self.locs))
            self.init_a.setall(False)
            self.init_a[:self.target] = True
        if self.init_a.count() != self.target:
            raise IndigoUnfeasibleComputation('Can only optimised when all '
                                'electrons are placed in the initial guess.')
        
    def run(self):
        self.initialise()
        min_ene = self.calc_energy(self.init_a)[1]
        seen = {self.init_a : min_ene}
        current_min = [self.init_a]
        min_round = min_ene + 1         # So the while loop is entered, 
        round_mins = current_min[:]     # regardless of min_ene value.
        pool = Pool(processes=NUM_PROCESSES)
        while abs(min_round - min_ene) > 1e-10:
            min_ene = min_round
            current_min = round_mins[:]
            a = current_min[0]
            results = pool.imap_unordered(self.calc_energy, 
                                               (n for n in self.neighbours(a)
                                                if n not in seen),
                                               chunksize=8)
            for n, n_ene in results:
                seen[n] = n_ene
                if n_ene - min_round < -1e-10:
                    min_round = n_ene
                    round_mins = [n]
                elif -1e-10 < n_ene - min_round < 1e-10:
                    round_mins.append(n)
        pool.terminate()
        bitarray_to_assignment(self.init_G, current_min[0], self.locs)
        return self.init_G, seen[current_min[0]]

    def calc_energy(self, a):
        graph_setup(self.G, a, self.locs)
                
        ene = sum(node_energy(self.G, n) for n in self.G)
        return a, round(ene, 5)
    
    def neighbours(self, a):
        for source in set(self.locs):
            i = self.locs.index(source)
            i_count = self.locs.count(source)
            i_loc = i + a[i:i+i_count].count() - 1
            if not a[i:i+i_count].count():
                continue
            for target in set(self.locs):
                if source == target:
                    continue
                j = self.locs.index(target)
                j_count = self.locs.count(target)
                j_loc = j + a[j:j+j_count].count()
                if j_count == a[j:j+j_count].count():
                    continue
                b = a.copy()
                b[i_loc] = False
                b[j_loc] = True
                yield b
                
    