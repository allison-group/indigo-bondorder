//
//  molecular_graph.cpp
//  indigox
//
//  Created by Welsh, Ivan on 25/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <iostream>
#include <map>
#include <tuple>

#include <boost/graph/connected_components.hpp>

#include "classes/atom.hpp"
#include "classes/bond.hpp"
#include "classes/iters.hpp"
#include "classes/molecular_graph.hpp"
#include "classes/molecule.hpp"
#include "utils/common.hpp"
#include "utils/options.hpp"

using namespace indigox;

MolecularGraph::MolecularGraph()
: _MolGraph()
{
}

MolecularGraph::MolecularGraph(Molecule_p m) : _MolGraph() {
  m->ResetIndices();
  std::map<Atom_p, MolVertex> atom_to_vertex;
  for (MolAtomIterator it = m->BeginAtom(); it != m->EndAtom(); ++it) {
    MolVertex v = AddVertex(*it);
    atom_to_vertex.emplace(*it, v);
  }
  
  for (MolBondIterator it = m->BeginBond(); it != m->EndBond(); ++it) {
    MolVertex u = atom_to_vertex.at((*it)->GetSourceAtom());
    MolVertex v = atom_to_vertex.at((*it)->GetTargetAtom());
    AddEdge(u, v, *it);
  }
}

MolVertex MolecularGraph::AddVertex(Atom_p atm) {
  MolVertProp p;
  p.atom = atm;
  return AddVertex(p);
}

MolEdgeBool MolecularGraph::AddEdge(MolVertex u, MolVertex v, Bond_p bnd) {
  MolEdgeProp p;
  p.bond = bnd;
  return AddEdge(u, v, p);
}

Uint MolecularGraph::NumConnectedComponents() {
  MolVertIdxMap idxMap;
  boost::associative_property_map<MolVertIdxMap> indexMap(idxMap);
  MolVertIterPair mvp = GetVertices();
  for (int i = 0; mvp.first != mvp.second; ++mvp.first, ++i) {
    boost::put(indexMap, (*mvp.first), i);
  }
  
  num_components_ = boost::connected_components(*graph_,
                                                boost::get(&MolVertProp::component, *graph_),
                                                boost::vertex_index_map(indexMap));
  return num_components_;
}

String MolecularGraph::ToDGFString() {
  typedef Options::AssignElectrons::FPT fpt_;
  std::ostringstream dgf_stream;
  uid_t u, v;
  MolEdgeIterPair es = GetEdges();
  for (auto start = es.first; start != es.second; ++start) {
    u = GetVertexIndex(GetSource(*start));
    v = GetVertexIndex(GetTarget(*start));
    if (fpt_::ADD_EDGES_TO_TD) {
      dgf_stream << "e " << u << "," << u << " " << u << "," << v << std::endl;
      dgf_stream << "e " << v << "," << v << " " << u << "," << v << std::endl;
    } else {
      dgf_stream << "e " << u << " " << v << std::endl;
    }
  }
  return dgf_stream.str();
}

