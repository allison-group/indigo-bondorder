//
//  electron_graph.hpp
//  indigox
//
//  Created by Welsh, Ivan on 12/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ELECTRON_GRAPH_HPP
#define ELECTRON_GRAPH_HPP

#include <cstdint>
#include <iostream>

#include "molecular_graph.hpp"
#include "../utils/graph.hpp"


namespace indigox {
  
  enum class SortOrder;
  
  struct ElnVertProp {
    MolVertPair id;
    float electronegativity;
    uint8_t valence;
    uint8_t atomic_number;
    uint8_t electron_count;
    uint8_t pre_placed;
    uint8_t target_octet;
    uint8_t target_hyper_octet;
    int8_t formal_charge;
    SortOrder sort_score;
  };
  
  typedef utils::Graph<ElnVertProp> _ElnGraph;
  typedef _ElnGraph::VertType ElnVertex;
  typedef _ElnGraph::VertIter ElnVertexIter;
  typedef _ElnGraph::EdgeType ElnEdge;
  typedef _ElnGraph::EdgeIter ElnEdgeIter;
  typedef _ElnGraph::NbrsIter ElnNeighboursIter;
  
  typedef _ElnGraph::VertIterPair ElnVertIterPair;
  typedef _ElnGraph::EdgeIterPair ElnEdgeIterPair;
  typedef _ElnGraph::NbrsIterPair ElnNbrsIterPair;
  
  
  class ElectronGraph : public _ElnGraph
  {
    
  public:
    ElectronGraph();// = default;
    ElectronGraph(const MolecularGraph &G);
    
    ElnVertex GetVertex(MolVertPair id) const;
  };
}

#endif /* ELECTRON_GRAPH_HPP */
