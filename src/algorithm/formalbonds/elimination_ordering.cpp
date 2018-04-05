//
//  elimination_ordering.cpp
//  indigox
//
//  Created by Ivan Welsh on 14/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//
#include <algorithm>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

#include "api.hpp"
#include "algorithm/formalbonds/elimination_ordering.hpp"
#include "classes/permutablegraph.hpp"

#ifdef BUILD_JAVA
#include "utils/java_interface.hpp"
#endif

namespace indigox {
  namespace algorithm {
    
    size_t VirtualElimination(PermutableGraph_p G, PermVertex v) {
      size_t addEdges = 0;
      PermNbrsIterPair outer = G->GetNeighbours(v);
      for (; outer.first != outer.second; ++outer.first) {
        PermNbrsIterPair inner = G->GetNeighbours(v);
        for (; inner.first != inner.second; ++inner.first) {
          if (*inner.first >= *outer.first) continue;
          if (!G->GetEdge(*inner.first, *outer.first).second) ++addEdges;
        }
      }
      return addEdges;
    }
    
    void RandomOrder(PermutableGraph_p G, ElimOrder& order) {
      order.clear();
      
      for (PermVertIterPair vs = G->GetVertices(); vs.first != vs.second; ++vs.first) {
        order.push_back(*vs.first);
      }
      
      std::random_device rd;
      std::mt19937 ran(rd());
      
      std::shuffle(order.begin(), order.end(), ran);
      
    }
    
    void QuickBBOrder(PermutableGraph_p G, ElimOrder& order) {
      order.clear();
      
#ifndef BUILD_JAVA
      std::cout << "To use the QuickBB algorithm, indigoX needs to be built with Java. Will fall back to using the Random algorithm. For shits and giggles." << std::endl;
      RandomOrder(G, order);
#else
      String DGFString = G->ToDGFString();
      String EO = utils::GetEliminationOrdering(DGFString);
      
      std::istringstream ss(EO);
      String idx;
      while (std::getline(ss, idx, ' ')) {
        indigox::uid_t i = (indigox::uid_t)std::stoull(idx);
        order.push_back(G->GetVertexByIndex(i).first);
      }
      
#endif
    }
    
    void MinDegreeOrder(PermutableGraph_p G, ElimOrder& order) {
      PermutableGraph_p g = PermutableGraph_p(new PermutableGraph(G));
      while (g->NumVertices() > 0) {
        size_t minD = g->NumVertices();
        PermVertex minV = NULL;
        PermVertIterPair verts = g->GetVertices();
        for (; verts.first != verts.second; ++verts.first) {
          if (g->Degree(*verts.first) < minD) {
            minD = g->Degree(*verts.first);
            minV = *verts.first;
          }
        }
        uid_t minIdx = g->GetVertexIndex(minV);
        order.push_back(G->GetVertexByIndex(minIdx).first);
        g->EliminateVertex(minV);
      }
    }
    
    void MinAddEdgesOrder(PermutableGraph_p G, ElimOrder& order) {
      PermutableGraph_p g = PermutableGraph_p(new PermutableGraph(G));
      while (g->NumVertices() > 0) {
        size_t minA = G->NumVertices();
        PermVertex minV = NULL;
        PermVertIterPair verts = g->GetVertices();
        for (; verts.first != verts.second; ++verts.first) {
          size_t currA = VirtualElimination(g, *verts.first);
          if (currA == 0) {
            minV = *verts.first;
          } else if (currA < minA) {
            minA = currA;
            minV = *verts.first;
          }
        }
        uid_t minIdx = g->GetVertexIndex(minV);
        order.push_back(G->GetVertexByIndex(minIdx).first);
        g->EliminateVertex(minV);
      }
    }
    
  }
}
