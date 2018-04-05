//
//  tree_decomposition.cpp
//  indigox
//
//  Created by Welsh, Ivan on 11/01/18.
//  Copyright © 2017 Allison Group. All rights reserved.
//
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>


#include "classes/molecular_graph.hpp"
#include "classes/permutablegraph.hpp"
#include "classes/treedecomp.hpp"

#include "utils/options.hpp"

#ifdef BUILD_JAVA
#include "utils/java_interface.hpp"
#endif

namespace indigox {
  
  
  typedef Options::AssignElectrons::FPT fpt_;
  
  TDecomp::TDecomp() : _TDGraph() { }
  
  TDecomp::TDecomp(PermutableGraph_p G, ElimOrder &order)
  : _TDGraph(), originalG_(G) {
    SetInput(G, order);
  }
  
  String TDecomp::ToString() {
    std::ostringstream ss;

    ss << "Bags: " << std::endl;
    TDVertIterPair bags = GetVertices();
    for (; bags.first != bags.second; ++bags.first) {
      ss << "    bag" << GetVertexIndex(*bags.first) << " : ";
      TDVertProp* p = GetProperties(*bags.first);
      for (PermVertex v : p->bag) {
        ss << originalG_->GetVertexIndex(v) << " ";
      }
      ss << std::endl;
    }
    ss << std::endl << "Edges: " << std::endl;
    TDEdgeIterPair edges = GetEdges();
    for (; edges.first != edges.second; ++edges.first) {
      TDVertex u = GetSource(*edges.first);
      TDVertex v = GetTarget(*edges.first);
      ss << "    bag" << GetVertexIndex(u) << "  --  ";
      ss << "bag" << GetVertexIndex(v) << std::endl;
    }
    return ss.str();
  }
  
  void TDecomp::SetInput(PermutableGraph_p G, ElimOrder &order) {
    Clear();
    originalG_ = PermutableGraph_p(G);
    permG_ = PermutableGraph_p(new PermutableGraph(G));
    uid_t eIdx = 0;
    for (PermVertex v : order) elimiIndex_.emplace(v, eIdx++);
    FromOrder(order, 0);
  }
  
  void TDecomp::FromOrder(ElimOrder &order, size_t index) {
    size_t s = permG_->NumVertices();
    
    if (s == 0) { // shouldn't happen
      upperBound_ = 0;
    } else if (s == 1) {
      upperBound_ = 0;
      TDVertProp p;
      p.bag.emplace(order[index]);
      AddVertex(p);
    } else if (s == 2) {
      upperBound_ = 1;
      TDVertProp p;
      p.bag.emplace(order[index]);
      p.bag.emplace(order[index + 1]);
      TDVertex v = AddVertex(p);
      PermVertIter it = permG_->GetVertices().first;
      bagMap_.emplace(*it++, v);
      bagMap_.emplace(*it++, v);
    } else {
      uid_t thisIdx = originalG_->GetVertexIndex(order[index]);
      PermVertex thisV = permG_->GetVertexByIndex(thisIdx).first;
      size_t numNbrs = permG_->Degree(thisV);
      
      uid_t lowIdx = std::numeric_limits<uid_t>::max();
      PermVertex lowNbr = NULL;
      TDVertProp p;
      p.bag.emplace(order[index]);
      PermNbrsIterPair nbrs = permG_->GetNeighbours(thisV);
      for (; nbrs.first != nbrs.second; ++nbrs.first) {
        uid_t currentIdx = permG_->GetVertexIndex(*nbrs.first);
        PermVertex originalV = originalG_->GetVertexByIndex(currentIdx).first;
        uid_t eIdx = elimiIndex_.at(originalV);
        p.bag.emplace(originalV);
        if (eIdx < lowIdx) {
          lowIdx = eIdx;
          lowNbr = *nbrs.first;
        }
      }
      
      permG_->EliminateVertex(thisV);
      FromOrder(order, index + 1);
      
      upperBound_ = std::max(upperBound_, numNbrs);
      TDVertex v = AddVertex(p);
      bagMap_.emplace(thisV, v);
      if (lowNbr != NULL) AddEdge(v, bagMap_.at(lowNbr));
    }
  }
}

