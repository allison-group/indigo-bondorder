//
//  molecular_graph.hpp
//  indigox
//
//  Created by Welsh, Ivan on 25/08/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP
#define INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP

#include <cstdint>
#include <map>
#include <string>
#include <tuple>

#include <boost/logic/tribool.hpp>

#include "../api.hpp"
#include "periodictable.hpp"
#include "../utils/counter.hpp"
#include "../utils/graph.hpp"

namespace indigox {
  
  struct MolVertProp {
//    float charge = 0.0f;
//    Element_p element;
//    int8_t formal_charge = 0;
//    // R,S as per normal, A for achiral, I for indetermined
//    char chirality = 'I';
//    bool aromaticity = false;
    Atom_p atom;
    int component = -1;
  };
  
  struct MolEdgeProp {
//    uint8_t bond_order = 1;
//    // E,Z as per normal, N for non-geometric bond, I for indetermined
//    char geometry = 'I';
    Bond_p bond;
  };
  
  // TODO: Replace OutEdgeList and VertexList with boost::vecS. Doing so
  //       will require handling of all invalidated descriptors on removal
  typedef utils::Graph<MolVertProp, MolEdgeProp> _MolGraph;
  typedef _MolGraph::VertType MolVertex;
  typedef _MolGraph::VertIter MolVertexIter;
  typedef _MolGraph::EdgeType MolEdge;
  typedef _MolGraph::EdgeIter MolEdgeIter;
  typedef _MolGraph::NbrsIter MolNeighboursIter;
  
  typedef _MolGraph::VertTypePair MolVertPair;
  typedef _MolGraph::VertIterPair MolVertIterPair;
  typedef _MolGraph::EdgeIterPair MolEdgeIterPair;
  typedef _MolGraph::NbrsIterPair MolNbrsIterPair;
  typedef _MolGraph::VertBool MolVertBool;
  typedef _MolGraph::EdgeBool MolEdgeBool;
  typedef std::tuple<MolVertex, MolVertex, MolEdgeProp> MolVertPairEdgeProp;
  
  typedef std::map<MolVertex, int> MolVertIdxMap;
  typedef std::map<MolEdge, int> MolEdgeIdxMap;
  
  class MolecularGraph : public _MolGraph {
  private:
    boost::tribool planar_ = boost::indeterminate;
    int16_t totalCharge_ = 0;
    Uint num_components_ = 1;
    Molecule_p source_;
    
  public:
    using _MolGraph::AddVertex;
    using _MolGraph::GetEdge;
    using _MolGraph::AddEdge;
    
    MolecularGraph();
    MolecularGraph(Molecule_p);
    
    MolVertex AddVertex(Atom_p);
    MolEdgeBool AddEdge(MolVertex, MolVertex, Bond_p);
    
    // Graph operations
    inline int16_t GetTotalCharge() const { return totalCharge_; }
    inline void SetTotalCharge(int16_t q) { totalCharge_ = q; }
    
    /// @todo transfer connected components to underlying graph
    Uint NumConnectedComponents();
    
    String ToDGFString();
    
  };
  
}  // namespace indigox

#endif /* INDIGOX_CLASSES_MOLECULAR_GRAPH_HPP */
