//
//  permutablegraph.cpp
//  indigox
//
//  Created by Ivan Welsh on 13/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//
#include <map>
#include <sstream>

#include "api.hpp"
#include "classes/molecular_graph.hpp"
#include "classes/permutablegraph.hpp"
#include "classes/treedecomp.hpp"
#include "utils/graph.hpp"
#include "utils/options.hpp"

using namespace indigox;
typedef Options::AssignElectrons::FPT fpt_;;

String PermutableGraph::ToDGFString() {
  std::ostringstream dgf;
  uid_t u, v;
  PermEdgeIterPair edges = GetEdges();
  while (edges.first != edges.second) {
    u = GetVertexIndex(GetSource(*edges.first));
    v = GetVertexIndex(GetTarget(*edges.first));
    dgf << "e " << u << " " << v << std::endl;
    ++edges.first;
  }
  return dgf.str();
}

String PermutableGraph::PGVToMGVTable() {
  std::ostringstream ss;
  PermVertIterPair verts = GetVertices();
  for (; verts.first != verts.second; ++verts.first) {
    ss << "PGV " << GetVertexIndex(*verts.first) << " --> MGV ";
    PermVertProp* p = GetProperties(*verts.first);
    uid_t u = source_->GetVertexIndex(p->source.first);
    uid_t v = source_->GetVertexIndex(p->source.second);
    if (u < v) ss << u << "," << v << std::endl;
    else ss << v << "," << u << std::endl;
  }
  return ss.str();
}

PermutableGraph::PermutableGraph() : _PermGraph() { }

PermutableGraph::PermutableGraph(MolecularGraph_p G) : _PermGraph() {
  SetInput(G);
}

PermutableGraph::PermutableGraph(PermutableGraph_p G) : _PermGraph() {
  source_ = G->GetSourceGraph();
  
  PermVertIterPair verts = G->GetVertices();
  std::map<PermVertex, PermVertex> originalToNew;
  while (verts.first != verts.second) {
    PermVertProp p;
    p.source = G->GetProperties(*verts.first)->source;
    PermVertex v = AddVertex(p);
    originalToNew.emplace(*verts.first, v);
    ++verts.first;
  }
  vert_idxmap_.clear();
  for (auto& o2n : originalToNew) {
    vert_idxmap_.insert({o2n.second, G->GetVertexIndex(o2n.first)});
  }
  
  PermEdgeIterPair edges = G->GetEdges();
  while (edges.first != edges.second) {
    PermVertex v1_p = G->GetSource(*edges.first);
    PermVertex v2_p = G->GetTarget(*edges.first);
    AddEdge(originalToNew.at(v1_p), originalToNew.at(v2_p));
    ++edges.first;
  }
}

void PermutableGraph::SetInput(MolecularGraph_p G) {
  Clear();
  source_ = G;
  
  MolEdgeIterPair edges = G->GetEdges();
  std::map<MolVertex, PermVertex> molToPerm;
  
  while (edges.first != edges.second) {
    MolVertex u = G->GetSource(*edges.first);
    MolVertex v = G->GetTarget(*edges.first);
    ++edges.first;
    
    PermVertProp p_u, p_v;
    PermVertex u_, v_;
    p_u.source = std::make_pair(u, u);
    p_v.source = std::make_pair(v, v);
    
    if (molToPerm.find(u) == molToPerm.end()) {
      u_ = AddVertex(p_u);
      molToPerm.emplace(u, u_);
    }
    else {
      u_ = molToPerm.at(u);
    }
    if (molToPerm.find(v) == molToPerm.end()) {
      v_ = AddVertex(p_v);
      molToPerm.emplace(v, v_);
    }
    else {
      v_ = molToPerm.at(v);
    }
    
    if (fpt_::ADD_EDGES_TO_TD) {
      PermVertProp p_w;
      if (u < v) p_w.source = std::make_pair(u, v);
      else p_w.source = std::make_pair(v, u);
      PermVertex w_ = AddVertex(p_w);
      AddEdge(u_, w_);
      AddEdge(v_, w_);
    } else {
      AddEdge(u_, v_);
    }
  }
}

void PermutableGraph::EliminateVertex(PermVertex v) {
  PermNbrsIterPair outer = GetNeighbours(v);
  while (outer.first != outer.second) {
    PermNbrsIterPair inner = GetNeighbours(v);
    while (inner.first != inner.second) {
      if (*inner.first != *outer.first) AddEdge(*inner.first, *outer.first);
      ++inner.first;
    }
    ++outer.first;
  }
  RemoveVertex(v);
}
