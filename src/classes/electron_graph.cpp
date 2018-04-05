//
//  electron_graph.cpp
//  indigox
//
//  Created by Welsh, Ivan on 12/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//
#include <iostream>
#include <unordered_map>

#include "classes/atom.hpp"
#include "classes/electron_graph.hpp"
#include "classes/molecular_graph.hpp"
#include "utils/graph.hpp"
#include "utils/options.hpp"

namespace indigox {
  
  uint8_t PreplaceCount(const MolecularGraph &G, MolVertex v) {
    if (!Options::AssignElectrons::PREPLACE_ELECTRONS) return 0;
    MolVertProp* p = G.GetProperties(v);
    switch (G.Degree(v)) {
      case 1:
        switch (p->atom->GetElement()->GetAtomicNumber()) {
          case 9:
          case 17:
          case 35:
            return 6;
          case 8:
          case 16:
            return 4;
          case 7:
            return 2;
          default:
            return 0;
        }
      case 2:
        switch (p->atom->GetElement()->GetAtomicNumber()) {
          case 8:
          case 16:
            return 4;
          default:
            return 0;
        }
      default:
        return 0;
    }
  }
  
  uint8_t PreplaceCount(const MolecularGraph &G, MolEdge e) {
    if (!Options::AssignElectrons::PREPLACE_ELECTRONS) return 2;
    MolVertex u = G.GetSource(e);
    MolVertex v = G.GetTarget(e);
    MolVertProp* u_p = G.GetProperties(u);
    MolVertProp* v_p = G.GetProperties(v);
    if (G.Degree(u) == 2 && G.Degree(v) == 1) {
      if (u_p->atom->GetElement() == 6 && v_p->atom->GetElement() == 7) {
        return 6;
      }
    } else if (G.Degree(v) == 2 && G.Degree(u) == 1) {
      if (u_p->atom->GetElement() == 7 && v_p->atom->GetElement() == 6) {
        return 6;
      }
    }
    return 2;
  }
  
  /************************
   *                      *
   *     Construction     *
   *                      *
   ************************/
  ElectronGraph::ElectronGraph()
  : _ElnGraph()
  {
  }
  
  ElectronGraph::ElectronGraph(const MolecularGraph &G)
  : _ElnGraph()
  {
    using namespace std;
    unordered_map<MolVertex, ElnVertex> mol2eln;
    MolVertIterPair vertices = G.GetVertices();
    uint16_t count = 0;
    for (MolVertexIter b = vertices.first; b != vertices.second; ++b) {
      MolVertProp prop = *G.GetProperties(*b);
      ElnVertProp e_prop;
      e_prop.id = make_pair(*b, *b);
      count++;
      e_prop.electronegativity = prop.atom->GetElement()->GetElectronegativity();
      e_prop.valence = prop.atom->GetElement()->GetValenceElectronCount();
      e_prop.atomic_number = prop.atom->GetElement()->GetAtomicNumber();
      e_prop.electron_count = 0;
      e_prop.pre_placed = PreplaceCount(G, *b);
      e_prop.target_octet = prop.atom->GetElement()->GetOctet();
      e_prop.target_hyper_octet = prop.atom->GetElement()->GetHypervalentOctet();
      e_prop.formal_charge = 0;
      
      
      ElnVertex e_vert = AddVertex(e_prop);
      mol2eln.insert(make_pair(*b, e_vert));
    }
    
    MolEdgeIterPair edges = G.GetEdges();
    for (MolEdgeIter b = edges.first; b != edges.second; ++b) {
      MolVertex u = G.GetSource(*b), v = G.GetTarget(*b);
      ElnVertex u_eln = mol2eln.at(u), v_eln = mol2eln.at(v);
      ElnVertProp e_prop;
      if (u < v)
        e_prop.id = make_pair(u, v);
      else
        e_prop.id = make_pair(v, u);
      e_prop.electron_count = 0;
      e_prop.pre_placed = PreplaceCount(G, *b);
      ElnVertex e_vert = AddVertex(e_prop);
      AddEdge(u_eln, e_vert);
      AddEdge(v_eln, e_vert);
    }
  }
  
  ElnVertex ElectronGraph::GetVertex(MolVertPair id) const {
    ElnVertIterPair vertices = GetVertices();
    for (ElnVertexIter it = vertices.first; it != vertices.second; ++it) {
      ElnVertProp* p = GetProperties(*it);
      if (p->id == id)
        return *it;
    }
    return ElnVertex();
  }
}

