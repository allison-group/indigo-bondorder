from collections import deque, defaultdict
import copy
from itertools import count as _count
from itertools import product
import os
from subprocess import TimeoutExpired

from indigoX.config import (JAVA_PATH, LIBTW_PATH, WORK_DIR, INFINITY, 
                            ELECTRON_PAIRS, ALLOW_HYPERVALENT, HYPERPENALTY, 
                            BASIS_LEVEL, RUN_QBND, COUNTERPOISE_CORRECTED, 
                            TD_TIMEOUT, ALLOW_FALLBACK, MAX_TREEWIDTH)
from indigoX.data import atom_enes, bond_enes, qbnd_enes
from indigoX.exception import IndigoExternalProgramError, IndigoUnfeasibleComputation
from indigoX.misc import (BondOrderAssignment, graph_to_dist_graph, electron_spots, electrons_to_add,
                          locs_sort, formal_charge, random_string)
from indigoX.periodictable import PeriodicTable as PT
import networkx as nx
import subprocess as sp


HALF_INF = INFINITY / 2
INF_SQUARED = INFINITY ** 2
BSSE = int(not COUNTERPOISE_CORRECTED)

class FPT(BondOrderAssignment):
    def __init__(self, G):
        self.init_G = G
        self.workdir = WORK_DIR / "BondOrderAssignment"
        self.workdir.mkdir(parents=True, exist_ok=True)
        self.max_vert = 0, None, 0
        self.best_root = 0
        self.current_root = 0
        
    def initialise(self):
        self.G = graph_to_dist_graph(self.init_G)
        self.target = electrons_to_add(self.init_G)
        self.locs = locs_sort(electron_spots(self.init_G), self.G)
        self.choices = []
        for n in self.locs:
            n_count = self.locs.count(n)
            if (n,n_count) not in self.choices:
                self.choices.append((n,n_count))
        for i in range(len(self.choices)):
            self.choices[i] = self.choices[i][1]

    def run(self):
        self.initialise()
        self.tree_decomposition()
        self.nice_tree_decomposition()
        self.setup_ntd()
        
        delete_after = defaultdict(list)
        for vert in self.vert_order:
            ntd_dat = self.NTD.node[vert]
            tmp = self.NTD.node[vert]
            tmp['score'] = TDVertScore(ntd_dat, self.G, self.locs, self.target)
            tmp['score'].propagate_with(*[self.NTD.node[n]['score'] 
                                          for n in self.NTD.successors(vert)])
            try:
                parent = self.NTD.predecessors(vert).pop()
            except IndexError:
                continue
            delete_after[parent].append(vert)
            for to_del in delete_after[vert]:
                del self.NTD.node[to_del]['score']
        
        f_dat = self.NTD.node[self.root]['score'].s_matrix[self.target]['final']
        self.forgotten_to_assignment(f_dat['forgotten'])   
        return (self.init_G, round(f_dat['energy'],5))
    
    def forgotten_to_assignment(self, forgotten):
        locs = []
        for l, c in forgotten:
            while c:
                locs.append(l)
                c -= 1
        
        self.G = graph_to_dist_graph(self.init_G)
        e_per_count = 2 if ELECTRON_PAIRS else 1
        for i in locs:
            self.G.node[i]['e-'] += e_per_count
        
        for n in self.G:
            if len(n) == 1:
                self.init_G.node[n[0]]['formal_charge'] = formal_charge(self.G, n)
            if len(n) == 2:
                self.init_G[n[0]][n[1]]['order'] = self.G.node[n]['e-'] // 2
            
    def setup_ntd(self):
        for parent in self.NTD:
            tmp = self.NTD.node[parent]
            if self.NTD.in_degree(parent) == 0:
                p_bag = tmp['bag']
                c_bag = self.NTD.node[self.NTD.neighbors(parent)[0]]['bag']
                tmp['kind'] = 'R', (c_bag - p_bag).pop()
            elif self.NTD.out_degree(parent) == 0:
                tmp['kind'] = 'L', None
            elif self.NTD.out_degree(parent) == 2:
                tmp['kind'] = 'J', None
            else:
                p_bag = tmp['bag']
                c_bag = self.NTD.node[self.NTD.neighbors(parent)[0]]['bag']
                if p_bag > c_bag:
                    tmp['kind'] = 'I', (p_bag - c_bag).pop()
                else:
                    tmp['kind'] = 'F', (c_bag - p_bag).pop()
                     
        paths = {}
        f_path = defaultdict(set)
        c = _count()
        for vert in self.NTD:
            if self.NTD.node[vert]['kind'][0] not in ['J','L']:
                continue
            j = next(c)
            path = nx.shortest_path(self.NTD, self.root, vert)
            path.reverse()
            for i in range(1,len(path)):
                if self.NTD.node[path[i]]['kind'][0] == 'J':
                    path = path[:i+1]
                    break
                elif self.NTD.node[path[i]]['kind'][0] in ['F','R']:
                    f_path[j].add(self.NTD.node[path[i]]['kind'][1])
             
            paths[j] = tuple(path)
         
        path_ordering = nx.DiGraph()    
        for k, path in paths.items():
            for k2, path2 in paths.items():
                if path2[-1] == path[0]:
                    path_ordering.add_edge(k, k2)
                 
        self.vert_order = []              
        for k in nx.topological_sort(path_ordering, reverse=True):
            if self.NTD.node[paths[k][-1]]['kind'][0] == 'R':
                self.vert_order += paths[k]
            else:
                self.vert_order += paths[k][:-1]
                 
        if set(self.vert_order) != set(self.NTD.nodes()):
            self.vert_order = nx.topological_sort(self.NTD, reverse=True)
                 
        for v in self.NTD:
            n_tot = self.target
            requires = nx.descendants(self.NTD, v).union({v})
            forgots = set(self.NTD.node[x]['kind'][1] for x in requires 
                          if self.NTD.node[x]['kind'][0] in ['F','R'])
            n_pos = sum(self.locs.count(x) for x in forgots)
            n_pos_2 = len(self.locs) - n_pos
            self.NTD.node[v]['min'] = max(0, n_tot - n_pos_2)
            self.NTD.node[v]['max'] = min(n_tot, n_pos)
                   
    def nice_tree_decomposition(self):
        c = _count(self.TD.order() + 1)
        # select a bag as root and make into directed tree
        root = self.possible_roots.pop()
        self.root = root
        self.NTD = nx.DiGraph()
        self.NTD.add_node(root, self.TD.node[root])
        p_c_verts = deque((root, child) for child in self.TD.neighbors(root))
        while p_c_verts:
            parent, child = p_c_verts.popleft()
            if child not in self.NTD:
                self.NTD.add_node(child, self.TD.node[child])
            self.NTD.add_edge(parent, child)
            for grandchild in self.TD.neighbors(child):
                if grandchild not in self.NTD:
                    p_c_verts.append((child, grandchild))
      
        # explicitly add edges into bags. Edges are introduced the vertex after
        # the first atom of the edge is added, and forgotten the vertex after
        # the second atom of the edge is forgotten. Only required if using 
        # charged bond energies
        for leaf in (x for x in self.NTD if not self.NTD.out_degree(x)):
            in_edges = []
            while self.NTD.in_degree(leaf):
                tmp = []
                for u, v in in_edges:
                    if ((u,) not in self.NTD.node[leaf]['bag']
                        and (v,) not in self.NTD.node[leaf]['bag']):
                        continue
                    tmp.append((u,v))
                in_edges = tmp
                for u in self.NTD.node[leaf]['bag']:
                    if len(u) != 1:
                        continue
                    for v in self.G.neighbors(u):
                        if v not in in_edges:
                            in_edges.append(v)
                     
                self.NTD.node[leaf]['bag'].update(in_edges)
                leaf = self.NTD.predecessors(leaf)[0]
 
        # reduce all vertices to only having 2 children
        while True:
            parent = None
            for v, deg in self.NTD.out_degree().items():
                if deg > 2:
                    parent = v
                    break
            else:
                break
            idx = next(c)
            self.NTD.add_node(idx, copy.deepcopy(self.NTD.node[parent]))
            child_to_remove = list(self.NTD.neighbors(parent))[1:]
            self.NTD.add_edge(parent, idx)
            for child in child_to_remove:
                self.NTD.remove_edge(parent, child)
                self.NTD.add_edge(idx, child)
                 
        # add same nodes after join nodes
        for parent in self.NTD.nodes():
            if self.NTD.out_degree(parent) != 2:
                continue
            child_a, child_b = self.NTD.neighbors(parent)
            if self.NTD.node[child_a]['bag'] != self.NTD.node[parent]['bag']:
                idx_a = next(c)
                self.NTD.add_node(idx_a, copy.deepcopy(self.NTD.node[parent]))
                self.NTD.add_edge(parent, idx_a)
                self.NTD.remove_edge(parent, child_a)
                self.NTD.add_edge(idx_a, child_a)
            if self.NTD.node[child_b]['bag'] != self.NTD.node[parent]['bag']:    
                idx_b = next(c)
                self.NTD.add_node(idx_b, copy.deepcopy(self.NTD.node[parent]))
                self.NTD.add_edge(parent, idx_b)
                self.NTD.remove_edge(parent, child_b)
                self.NTD.add_edge(idx_b, child_b)
         
        # add an intersection node between each parent/child where the child is
        # not a subset of the parent and the parent is not a subset of the child
        for parent in self.NTD.nodes():
            p_bag = self.NTD.node[parent]['bag']
            for child in self.NTD.neighbors(parent):
                c_bag = self.NTD.node[child]['bag']
                if p_bag.issubset(c_bag) or c_bag.issubset(p_bag):
                    continue
                i_bag = p_bag.intersection(c_bag)
                i_idx = next(c)
                self.NTD.add_node(i_idx, {'bag':i_bag})
                self.NTD.remove_edge(parent, child)
                self.NTD.add_edge(parent, i_idx)
                self.NTD.add_edge(i_idx, child)
         
        # add introduce and forget nodes as needed
        for parent in nx.topological_sort(self.NTD, reverse=True):
            p_bag = self.NTD.node[parent]['bag']
            for child in self.NTD.neighbors(parent):
                c_bag = self.NTD.node[child]['bag']
                if p_bag == c_bag:
                    continue
                elif p_bag.issuperset(c_bag):
                    while len(p_bag) > len(c_bag) + 1:
                        # introduce vertices with smallest loc count first
                        i_bag = sorted(p_bag.difference(c_bag), 
                                       key=lambda x:self.locs.count(x), 
                                       reverse=True).pop()
                        i_bag = c_bag.union({i_bag})
                        i_idx = next(c)
                        self.NTD.add_node(i_idx, {'bag':i_bag})
                        self.NTD.remove_edge(parent, child)
                        self.NTD.add_edge(parent, i_idx)
                        self.NTD.add_edge(i_idx, child)
                        c_bag = i_bag
                        child = i_idx
                elif p_bag.issubset(c_bag):
                    while len(p_bag) < len(c_bag) - 1:
                        # forget vertices with largest loc count first
                        f_bag = sorted(c_bag.difference(p_bag), 
                                       key=len, #lambda x:self.locs.count(x), 
                                       reverse=True).pop()
                        f_bag = c_bag.difference({f_bag})
                        f_idx = next(c)
                        self.NTD.add_node(f_idx, {'bag':f_bag})
                        self.NTD.remove_edge(parent, child)
                        self.NTD.add_edge(parent, f_idx)
                        self.NTD.add_edge(f_idx, child)
                        c_bag = f_bag
                        child = f_idx

    def tree_decomposition(self):
        # write the input file
        tmp_str = random_string(8, ('ABCDEFGHIJKLMNOPQRSTUVWXYZqwertyuioplkj'
                                    'hgfdsazxcvbnm1234567890'))
        in_file = self.workdir/(tmp_str+'.dgf')
        print('c',self.init_G.graph['name'], file=in_file.open('w'))
        for u, v in self.init_G.edges():
            edge_l = ['e']
            edge_l.append('{0},{0}'.format(u))
            edge_l.append('{0},{0}'.format(v))
            print(' '.join(edge_l), file=in_file.open('a'))
 
        td_args = [str(JAVA_PATH), '-jar', 'libtw.jar', 
                   'nl.uu.cs.treewidth.TDPrint', 'QuickBB', str(in_file)]
        TD_process = sp.Popen(td_args, stdout=sp.PIPE, cwd=str(LIBTW_PATH))
        try:
            td_str, _ = TD_process.communicate(timeout=TD_TIMEOUT)
        except TimeoutExpired:
            TD_process.kill()
            TD_process.communicate()
            if not ALLOW_FALLBACK:
                os.remove(str(in_file))
                raise IndigoExternalProgramError("QuickBB timeout ({} seconds)"
                                                 "".format(TD_TIMEOUT))
            # If timed out on QuickBB, run with upperbound only
            self.log.warning('QuickBB timeout ({} seconds). Falling back to '
                             'using GreedyFillIn upperbound as TD.'
                             ''.format(TD_TIMEOUT))
            td_args[-2] = 'upperbound'
            TD_process = sp.Popen(td_args, stdout=sp.PIPE, cwd=str(LIBTW_PATH))
            td_str, _ = TD_process.communicate()
            
        in_v = in_e = done_v = False
        self.TD = nx.Graph()
        self.tree_width = 0
        decoded_output = td_str.decode()
        first_line = decoded_output.split('\n')[0]
        if 'error' in first_line:
            os.remove(str(in_file))
            raise IndigoExternalProgramError("libtw.jar failed to provide a "
                                             "tree-decomposition with: {}"
                                             "".format(first_line))
            
        for line in decoded_output.split('\n'):
            line = line.strip()
            if line.startswith('bag') and not in_v:
                in_v = True
            elif line.startswith('bag') and in_v and not done_v:
                pass
            elif line.startswith('bag') and done_v and not in_e:
                in_e = True
            elif line.startswith('bag') and in_e:
                pass
            elif in_v and not done_v and not line:
                done_v = True
                continue
            elif in_e and not line:
                break
            
            if in_v and not done_v:
                idx = int(line.split()[0][3:])
                bag = set()
                for v in line.split('"')[1].split():
                    a, b = map(int, v.split(','))
                    if a == b:
                        bag.add((a,))
                    else:
                        bag.add((a,b))
                self.TD.add_node(idx, {'bag':bag})
                if len(bag) - 1 > self.tree_width:
                    self.tree_width = len(bag) - 1
            elif in_e:
                a_idx = int(line.split()[0][3:])
                b_idx = int(line.split()[2][3:])
                self.TD.add_edge(a_idx, b_idx)
                
        c = _count(self.TD.order() + 1)
        self.possible_roots = []
        # Add an empty bag to every leaf
        for n in self.TD.nodes():
            if self.TD.degree(n) == 1 or self.TD.degree(n) == 0:
                idx = next(c)
                self.TD.add_node(idx, {'bag':set()})
                self.TD.add_edge(n, idx)
                self.possible_roots.append(idx)
        
        # remove the temporary input file for TD
        os.remove(str(in_file))
        if MAX_TREEWIDTH > 0 and self.tree_width > MAX_TREEWIDTH:
            raise IndigoUnfeasibleComputation('Treewidth is too large')
        
class TDVertScore(object):
    def __init__(self, NTD_data, G, locs, global_target):
        self.bag = NTD_data['bag']
        self.locs = locs
        self.G = G
        self.s_matrix = defaultdict(lambda: defaultdict(dict))
        self.kind = NTD_data['kind'][0]
        self.change_v = NTD_data['kind'][1]
        self.global_max = global_target
        self.min_count = NTD_data['min']
        self.max_count = NTD_data['max']
        
        tmp = []
        for vert in sorted(self.bag):
            tmp_2 = []
            step = 2 if len(vert) == 2 and not ELECTRON_PAIRS else 1
            for i in range(0,locs.count(vert) + 1, step):
                tmp_2.append((vert,i))
            tmp.append(tmp_2)
        
        for forget_e in range(self.min_count, self.max_count + 1):
            for row in product(*tmp):
                if not row:
                    continue
                fs = frozenset(row)
                if fs not in self.s_matrix[forget_e]:
                    ene = 0 if self.kind not in ['J','F'] else INF_SQUARED
                    self.s_matrix[forget_e][fs] = {'forgotten':[], 
                                                   'energy':ene,
                                                   'others':[],
                                                   }
        if self.kind == 'R':
            self.s_matrix[self.global_max]['final'] = {'forgotten':[],
                                                       'energy':INF_SQUARED,
                                                       'others':[],
                                                       }
        
    
    def propagate_with(self, a=None, b=None):
        if self.kind == 'L':
            self._leaf(a, b)
        elif self.kind == 'R':
            self._root(a, b)
        elif self.kind == 'J':
            self._join(a, b)
        elif self.kind == 'I':
            self._introduce(a, b)
        elif self.kind == 'F':
            self._forget(a, b)
            
            
    def _leaf(self, a, b):
        # Nothing has to be done for the leafs as they're already covered
        # by the init
        pass
    
    def _root(self, a, b):
        # same as a forget?
        # Need to score the forgotten vertex for each row of s_matrix        
        for a_e in a.s_matrix:
            vert_scores = {}
            all_enes = {k:a.s_matrix[a_e][k]['energy'] for k in a.s_matrix[a_e]}
            all_forgets = {k:a.s_matrix[a_e][k]['forgotten'] for k in a.s_matrix[a_e]}
            for fs in a.s_matrix[a_e]:
                not_in_parent = None
                for v in fs:
                    if v[0] == self.change_v:
                        not_in_parent = v
                        break
                e = a_e + v[1]
                if e not in self.s_matrix:
                    continue
                if all_enes[fs] > HALF_INF:
                    continue
                if fs not in vert_scores:
                    vert_scores[fs] = self._score_vert(not_in_parent, fs.union(a.s_matrix[a_e][fs]['forgotten']))
                for parent_fs in self.s_matrix[e]:
                    forget = set(fs).pop()
                    ene = all_enes[fs] + vert_scores[fs]
                    if ene > HALF_INF:
                        ene = INFINITY
                    if self.s_matrix[e][parent_fs]['energy'] > ene:
                        self.s_matrix[e][parent_fs]['energy'] = ene
                        self.s_matrix[e][parent_fs]['forgotten'] = all_forgets[fs] + [forget]
                              
    def _join(self, a, b):
        # Need to combine the two forgottens and scores
        #print('------')
        for a_e in a.s_matrix:
            for b_e in b.s_matrix:
                e = a_e + b_e
                if e not in range(self.min_count, self.max_count+1):
                    continue
                for fs in self.s_matrix[e]:
                    comb_ene = a.s_matrix[a_e][fs]['energy'] + b.s_matrix[b_e][fs]['energy']
                    if comb_ene > HALF_INF:
                        comb_ene = INFINITY
                    if comb_ene < self.s_matrix[e][fs]['energy']:
                        self.s_matrix[e][fs]['energy'] = comb_ene
                        self.s_matrix[e][fs]['forgotten'] = a.s_matrix[a_e][fs]['forgotten'] + b.s_matrix[b_e][fs]['forgotten']
    
    def _introduce(self, a, b):
        # Need to propagate up the forgottens and energies
        
        for e in self.s_matrix:
            for fs in self.s_matrix[e]:
                not_in_child = None
                for v in fs:
                    if v[0] == self.change_v:
                        not_in_child = v
                        break
                child_fs = fs - {not_in_child}
                if child_fs in a.s_matrix[e]:
                    tmp = a.s_matrix[e][child_fs]
                    self.s_matrix[e][fs] = {'energy':tmp['energy'],
                                            'forgotten':tmp['forgotten'],
                                            'others':tmp['others'],
                                            }
                                
    def _score_vert(self, forgotten, others):
        v = forgotten[0]
        ene = INFINITY
        e_per_step = 2 if ELECTRON_PAIRS else 1
        if sum(x[1] for x in others.union({forgotten})) > self.global_max:
            return ene
        
        if len(v) == 1:
            element = self.G.node[v]['Z']
            octet = (PT[element].hyper if ALLOW_HYPERVALENT 
                     and self.G.degree(v) > 2 else PT[element].octet)
            fc = self.G.node[v]['valence'] - self.G.node[v]['e-']  # pre placed 
            fc -= forgotten[1] * e_per_step                        # me placed
            fc -= self.G.degree(v)                                 # bonded
            nbs = self.G.neighbors(v)
            val = self.G.node[v]['e-'] + forgotten[1] * e_per_step 
            val += 2 * self.G.degree(v)
            for n, e in others:
                if n not in nbs:
                    continue
                fc -= e // 2 if not ELECTRON_PAIRS else e
                val += e_per_step * e
                
            try:
                ene = atom_enes[BASIS_LEVEL][element][fc]
            except KeyError:
                pass
            
            if HYPERPENALTY and val > octet:
                ene = INFINITY
        elif len(v) == 2 and not RUN_QBND:
            a, b = v
            order = self.G.node[v]['e-'] + e_per_step * forgotten[1]
            a_dat = self.G.node[(a,)]
            b_dat = self.G.node[(b,)]
            a = a_dat['Z'] 
            b = b_dat['Z']
            a, b = sorted((a, b))
            if order % 2:
                ene = order / 2
            else:
                ene = bond_enes[BASIS_LEVEL][(a, b, order//2)][BSSE]
                
        elif len(v) == 2 and RUN_QBND:
            a, b = v
            order = self.G.node[v]['e-'] + e_per_step * forgotten[1]
            a_dat = self.G.node[(a,)]
            b_dat = self.G.node[(b,)]
            # Need to make sure don't violate valency of either of the atoms
            # in the bond
            a_fc = a_dat['valence'] - a_dat['e-'] - self.G.degree((a,))
            b_fc = b_dat['valence'] - b_dat['e-'] - self.G.degree((b,))
            a_nbs = self.G.neighbors((a,))
            b_nbs = self.G.neighbors((b,))
            for n, e in others:
                if n == (a,):
                    a_fc -= e * e_per_step
                elif n == (b,):
                    b_fc -= e * e_per_step
                if n in a_nbs:
                    a_fc -= e // 2 if not ELECTRON_PAIRS else e
                if n in b_nbs:
                    b_fc -= e // 2 if not ELECTRON_PAIRS else e
                
            a_sign = '+' if int(a_fc) > 0 else ''
            a_sign = '-' if int(a_fc) < 0 else a_sign
            b_sign = '+' if int(b_fc) > 0 else ''
            b_sign = '-' if int(b_fc) < 0 else b_sign
            a = a_dat['Z'] + a_sign 
            b = b_dat['Z'] + b_sign
            a, b = sorted([a,b])
            if order % 2:
                ene = order / 2
            elif (a, b, order//2) in qbnd_enes[BASIS_LEVEL]:
                ene = qbnd_enes[BASIS_LEVEL][(a, b, order//2)][BSSE]
            else:
                a = a[:-1] if ('+' in a or '-' in a) else a
                b = b[:-1] if ('+' in b or '-' in b) else b
                a,b = sorted((a,b))
                try:
                    ene = bond_enes[BASIS_LEVEL][(a, b, order//2)][BSSE]
                except KeyError:
                    pass    
        return ene
            
    def _forget(self, a, b):
        # Need to score the forgotten vertex for each row of s_matrix
        # for each forget electron count
        for a_e in a.s_matrix:
            vert_scores = {}
            all_enes = {k:a.s_matrix[a_e][k]['energy'] for k in a.s_matrix[a_e]}
            all_forgets = {k:a.s_matrix[a_e][k]['forgotten'] for k in a.s_matrix[a_e]}
            for fs in a.s_matrix[a_e]:
                not_in_parent = None
                for v in fs:
                    if v[0] == self.change_v:
                        not_in_parent = v
                        break
                e = a_e + v[1]
                if e not in self.s_matrix:
                    continue
                for parent_fs in self.s_matrix[e]:
                    if not fs.issuperset(parent_fs):
                        continue
                    if all_enes[fs] > HALF_INF:
                        continue
                    if fs not in vert_scores:
                        vert_scores[fs] = self._score_vert(not_in_parent, fs.union(a.s_matrix[a_e][fs]['forgotten']))
                    if vert_scores[fs] > HALF_INF:
                        continue
                    forget = set(fs - parent_fs).pop()
                    ene = all_enes[fs] + vert_scores[fs]
                    if ene > HALF_INF:
                        continue
                    if self.s_matrix[e][parent_fs]['energy'] > ene:
                        self.s_matrix[e][parent_fs]['energy'] = ene
                        self.s_matrix[e][parent_fs]['forgotten'] = all_forgets[fs] + [forget]
                        
