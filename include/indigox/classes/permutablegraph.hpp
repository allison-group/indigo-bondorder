//
//  permutablegraph.hpp
//  indigox
//
//  Created by Ivan Welsh on 13/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//

#ifndef INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP
#define INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP

#include <memory>
#include <vector>

#include "../api.hpp"
#include "molecular_graph.hpp"
#include "../utils/graph.hpp"

namespace indigox {
  struct PermVertProp {
    MolVertPair source;
    uid_t bag;
  };
  
  typedef utils::Graph<PermVertProp> _PermGraph;
  typedef _PermGraph::VertType PermVertex;
  typedef _PermGraph::NbrsIter PermNbrsIter;
  typedef _PermGraph::VertIter PermVertIter;
  typedef _PermGraph::NbrsIterPair PermNbrsIterPair;
  typedef _PermGraph::VertIterPair PermVertIterPair;
  typedef _PermGraph::EdgeIterPair PermEdgeIterPair;
  
  typedef std::vector<PermVertex> ElimOrder;
  
  class PermutableGraph : public _PermGraph
  {
  public:
    PermutableGraph();
    PermutableGraph(MolecularGraph_p);
    PermutableGraph(PermutableGraph_p);
    
  public:
    void SetInput(MolecularGraph_p);
    void EliminateVertex(PermVertex);
    String ToDGFString();
    MolecularGraph_p GetSourceGraph() { return source_; }
    String PGVToMGVTable();
    
  private:
    MolecularGraph_p source_;
  };
  
  typedef std::shared_ptr<PermutableGraph> PermutableGraph_p;
  
}

#endif /* INDIGOX_CLASSES_PERMUTABLEGRAPH_HPP */
