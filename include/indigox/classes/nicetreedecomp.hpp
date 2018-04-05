//
//  tree_decomposition.hpp
//  indigox
//
//  Created by Welsh, Ivan on 20/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_CLASSES_NICETREEDECOMP_HPP
#define INDIGOX_CLASSES_NICETREEDECOMP_HPP

#include <set>
#include <string>

#include "../api.hpp"
#include "../utils/graph.hpp"
#include "treedecomp.hpp"

namespace indigox {
    
    struct NTDVertProp {
      uint32_t id;
      std::set<MolVertPair> bag;
      std::pair<char, MolVertPair> kind;
    };
  
    
  typedef utils::Graph<NTDVertProp, utils::NoProperty,
  utils::DirectedGraph> _NTDGraph;
    
    typedef _NTDGraph::VertType NTDVertex;
    typedef _NTDGraph::VertIter NTDVertexIter;
    typedef _NTDGraph::EdgeType NTDEdge;
    typedef _NTDGraph::EdgeIter NTDEdgeIter;
    typedef _NTDGraph::NbrsIter NTDNeighboursIter;
    typedef _NTDGraph::PredIter NTDPredecessorsIter;
    
    typedef _NTDGraph::VertTypePair NTDVertPair;
    typedef _NTDGraph::VertIterPair NTDVertIterPair;
    typedef _NTDGraph::EdgeIterPair NTDEdgeIterPair;
    typedef _NTDGraph::NbrsIterPair NTDNbrsIterPair;
    typedef _NTDGraph::PredIterPair NTDPredIterPair;
    typedef _NTDGraph::VertBool NTDVertBool;
    typedef _NTDGraph::EdgeBool NTDEdgeBool;
    
    
    
    class NTDecomp : public _NTDGraph
    {
    public:
      NTDecomp();
      NTDecomp(TDecomp_p);
      NTDecomp(TDecomp_p, TDVertex);
      
      void SetInput(TDecomp_p);
      void SetInput(TDecomp_p, TDVertex);
      void TopologicalSort(std::vector<NTDVertex>&) const;
      TDecomp_p GetSourceGraph() { return source_; }
      size_t GetWidth();
      
    private:
      uint32_t vert_count_ = 0;
      TDecomp_p source_;
      
    };
  
}

#endif /* INDIGOX_CLASSES_NICETREEDECOMP_HPP */
