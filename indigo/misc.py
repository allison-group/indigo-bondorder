from abc import ABCMeta, abstractmethod
from bitarray import bitarray
from indigo.data import atom_enes, bond_enes, qbnd_enes, lp_prob, bo_prob
from indigo.config import (INFINITY, RUN_QBND, BASIS_LEVEL, HYPERVALENT, ELECTRON_PAIRS, HYPERPENALTY, 
                           PREFILL_LOCATIONS, COUNTERPOISE_CORRECTED)
from indigo.periodictable import PeriodicTable as PT
from indigo.exception import IndigoUnfeasibleComputation
from random import randrange
import networkx as nx
import logging.handlers
import sys

BSSE = int(not COUNTERPOISE_CORRECTED)

def node_energy(G, n):
    # Calculate the energy of a node in a BOAssignment graph
    ene = INFINITY
    
    if len(n) == 1:
        e = G.node[n]['Z']
        fc = G.node[n]['fc']
        es = G.node[n]['e-']
        try:
            ene = atom_enes[BASIS_LEVEL][e][fc]
        except KeyError:
            pass
        if HYPERPENALTY:
            val = es + sum(G.node[x]['e-'] for x in G.neighbors(n))
            octet = PT[e].hyper if HYPERVALENT and G.degree(n) > 2 else PT[e].octet
            if val > octet:
                ene = INFINITY

    elif len(n) == 2:
        a, b = n
        order = G.node[n]['e-']
        a_sign = '+' if int(G.node[(a,)]['fc']) > 0 else ''
        a_sign = '-' if int(G.node[(a,)]['fc']) < 0 else a_sign
        b_sign = '+' if int(G.node[(b,)]['fc']) > 0 else ''
        b_sign = '-' if int(G.node[(b,)]['fc']) < 0 else b_sign
        a = G.node[(a,)]['Z'] + a_sign 
        b = G.node[(b,)]['Z'] + b_sign
        a, b = sorted((a, b))
        if order % 2:
            ene = order / 2
        elif RUN_QBND and (a, b, order//2) in qbnd_enes[BASIS_LEVEL]:
            ene = qbnd_enes[BASIS_LEVEL][(a, b, order//2)][BSSE]
        else:
            a = a[:-1] if ('+' in a or '-' in a) else a
            b = b[:-1] if ('+' in b or '-' in b) else b
            try:
                ene = bond_enes[BASIS_LEVEL][(a, b, order//2)][BSSE]
            except KeyError:
                pass
            
    return ene

def formal_charge(G, a):
    # Calculates the formal charge on an atom given the bonding environment
    fc = G.node[a]['valence'] - G.node[a]['e-']
    for n in G.neighbors(a):
        e = G.node[n]['e-']
        za = G.node[a]['Z']
        zb = G.node[(n[0],) if n[1] == a[0] else (n[1],)]['Z']
        if not e % 2:
            fc -= e // 2
        # split electrons from odd count bonds when elements are the same
        elif e % 2 and za == zb and a[0] == n[0]:
            fc -= (e + 1) // 2
        elif e % 2 and za == zb:
            fc -= (e - 1) // 2
        # give odd electron from bond to most electronegative atom when
        # element different
        elif e % 2 and PT[za].chi > PT[zb].chi:
            fc -= (e + 1) // 2
        else:
            fc -= (e - 1) // 2
     
    return fc

def graph_to_dist_graph(G):
    # Convert the molecular graph to bond order assignment graph
    H = nx.Graph()
    for atom, dat in G.nodes(True):
        H.add_node((atom,), {'Z': dat['element'],
                             'e-': 0,
                             'valence': PT[dat['element']].valence,
                             'fc': 0,})

    for a, b, dat in G.edges(data=True):
        a, b = sorted((b, a))
        H.add_node((a, b), {'e-': 2})

        H.add_edge((a,), (a, b))
        H.add_edge((a, b), (b,))
        
    if PREFILL_LOCATIONS:
        for v, d in H.nodes(True):
            if len(v) == 2:
                continue
            H.node[v]['prefill'] = 0
            if lp_prob[(d['Z'], H.degree(v))][4] > 1000:
                for i in lp_prob[(d['Z'], H.degree(v))][1:4]:
                    if i == 1.0:
                        H.node[v]['prefill'] += 2 
    return H

def electron_spots(G):
    # Determines the places where electrons can be placed.
    spots = []
    for n, d in G.nodes(True):
        if HYPERVALENT and G.degree(n) > 2:
            octet = PT[d['element']].hyper
        else:
            octet = PT[d['element']].octet
            
        bonded = G.degree(n) * 2
        missing_e = octet - bonded
        if PREFILL_LOCATIONS:
            if lp_prob[(d['element'], G.degree(n))][4] > 1000:
                for i in lp_prob[(d['element'], G.degree(n))][1:4]:
                    if i == 1.0:
                        missing_e -= 2
            
        while missing_e > 0:
            spots.append((n,))
            if ELECTRON_PAIRS:
                missing_e -= 2
            else:
                missing_e -= 1
    for a, b in G.edges():
        a, b = sorted((a, b))
        if HYPERVALENT and G.degree(a) > 2:
            a_octet = PT[G.node[a]['element']].hyper
        else:
            a_octet = PT[G.node[a]['element']].octet
        if HYPERVALENT and G.degree(b) > 2:
            b_octet = PT[G.node[b]['element']].hyper
        else:
            b_octet = PT[G.node[b]['element']].octet
        a_bonded = G.degree(a) * 2
        b_bonded = G.degree(b) * 2
        a_missing = a_octet - a_bonded
        b_missing = b_octet - b_bonded
        order = 1
        a_e = G.node[a]['element']
        b_e = G.node[b]['element']
        a_e, b_e = sorted((a_e, b_e))
        while (a_missing > 0 and b_missing > 0
                and (a_e, b_e, order + 1) in bond_enes[BASIS_LEVEL]):
            spots.append((a, b))
            if not ELECTRON_PAIRS:
                spots.append((a, b))
            a_missing -= 2
            b_missing -= 2
            order += 1

    return spots

def electrons_to_add(G):
    # Determines how many electrons/electron pairs to add to the system
    total_e = sum(PT[G.node[n]['element']].valence for n in G)
    total_e -= G.graph['total_charge']
    total_e -= G.size() * 2         # Implicitly have all bonds with an electron
                                    # pair in them.
    if PREFILL_LOCATIONS:
        for n, d in G.nodes(True):
            if lp_prob[(d['element'], G.degree(n))][4] > 1000:
                for i in lp_prob[(d['element'], G.degree(n))][1:4]:
                    if i == 1.0:
                        total_e -= 2                           
    
    if ELECTRON_PAIRS and total_e % 2:
        raise IndigoUnfeasibleComputation('Unable to handle odd number of '
                                        'electrons when using electron pairs.')
    elif ELECTRON_PAIRS:
        total_e = total_e // 2
        
    return total_e 

def locs_sort(locs, G):
    # sorts the possible locations based on the probability that they will be
    # filled. Most probable first.
    def probability_sort(n, G=G):
        if len(n) == 1:
            e = G.node[n]['Z']
            d = G.degree(n)
            return lp_prob[(e,d)]
        elif len(n) == 2:
            a_e = G.node[(n[0],)]['Z']
            b_e = G.node[(n[1],)]['Z']
            a_d = G.degree((n[0],))
            b_d = G.degree((n[1],))
            if a_e > b_e:
                a_e, b_e = b_e, a_e
                a_d, b_d = b_d, a_d
            elif a_e == b_e and b_d < a_d:
                a_d, b_d = b_d, a_d
            return bo_prob[(a_e, a_d, b_e, b_d)]
        
    everything = set(locs)
    
    everything = sorted(everything, key=probability_sort, reverse=True)
    new_order = []
    
    for n in everything:
        while n in locs:
            locs.remove(n)
            new_order.append(n) 

    return new_order

def graph_setup(G, a, locs):
    # Setsup the BOAssign graph for heuristic and actual energy calculation
    for n in G:
        if len(n) == 1 and not PREFILL_LOCATIONS:
            G.node[n]['e-'] = 0
        elif len(n) == 1 and PREFILL_LOCATIONS:
            G.node[n]['e-'] = G.node[n]['prefill']
        elif len(n) == 2:
            G.node[n]['e-'] = 2
            
    for i in range(a.length()):
        if ELECTRON_PAIRS and a[i]:
            G.node[locs[i]]['e-'] += 2
        elif a[i]:
            G.node[locs[i]]['e-'] += 1
            
    for n in G:
        if len(n) == 1:
            G.node[n]['fc'] = formal_charge(G, n)
            
def bitarray_to_reallocs(a, locs):
    # Converts a bitarray into the locations it corresponds to
    r_locs = []
    for i in range(a.length()):
        if a[i]:
            r_locs.append(locs[i])
    return sorted(r_locs)

def calculable_nodes(G, a, stop, locs, target):
    # Determines which nodes of a BOAssign graph are calculable
    more_locs = locs[stop:]
    if a.count() >= target:
        return set(G.nodes())

    calculable = []
    for n in G:
        G.node[n]['changeable'] = True if n in more_locs else False
        
    for n in sorted(G, key=len):
        if G.node[n]['changeable']:
            continue
        if len(n) == 1:
            # atoms are calculable if they and their neighbours are unchangeable
            for nb in G[n]:
                if G.node[nb]['changeable']:
                    break
            else:
                calculable.append(n)
                
        elif len(n) == 2:
            # bonds are calculable if they and their neighbours are unchangeable
            for nb in G[n]:
                if nb not in calculable:
                    break
                elif G.node[nb]['changeable']:
                    break
            else:
                calculable.append(n)
      
    return set(calculable)

def random_string(length=4, chars="ABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"):
    return ''.join(chars[randrange(0,len(chars))] for _ in range(length))

def _get_logger():
    VERBOSE = 5 
    logging.addLevelName(VERBOSE, 'VERBOSE')
    logging.Logger.verbose = (lambda inst, msg, *args, **kwargs: 
                          inst.log(VERBOSE, msg, *args, **kwargs))
    logging.captureWarnings(True)
    _simple = logging.Formatter("[{levelname:5}]: {message}",style='{')
    _with_time = logging.Formatter("[{levelname:5}] {threadName}.{processName} "
                              "{module}.{funcName}.{lineno}: {message}\n"
                              "Time: {relativeCreated:.6f} ms", 
                              style='{')
    _long = logging.Formatter("[{levelname:5}] {threadName}.{processName} "
                              "{module}.{funcName}.{lineno}: {message}", 
                              style='{')

    _log_level = 20
    _c_handler = logging.StreamHandler(stream=sys.stderr)
    _c_handler.setLevel(_log_level)
    _c_handler.setFormatter(_simple)
    log = logging.getLogger()
    log.setLevel(_log_level)
    log.addHandler(_c_handler)
    return log 

class HashBitArray(bitarray):
    # A hashable version of bitarray.
    def __hash__(self):
        return hash(self.to01())
    
    def __eq__(self, x):
        if not isinstance(x, bitarray):
            return False
        elif self.length() != x.length():
            return False
        else:
            return self.to01() == x.to01()

class BondOrderAssignment(object):
    __metaclass__ = ABCMeta
    log = _get_logger()
    
    @abstractmethod
    def __init__(self, G, *args):
        pass
    
    @abstractmethod
    def run(self):
        pass
    
    @abstractmethod
    def initialise(self):
        pass
