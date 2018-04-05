//
//  tree_decomposition.cpp
//  indigox
//
//  Created by Welsh, Ivan on 20/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//
#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <boost/graph/topological_sort.hpp>

#include "api.hpp"

#include "classes/molecular_graph.hpp"
#include "classes/nicetreedecomp.hpp"
#include "classes/treedecomp.hpp"
#include "utils/options.hpp"

namespace indigox {
  
  size_t NTDecomp::GetWidth() {
    size_t maxBag = 0;
    NTDVertIterPair verts = GetVertices();
    for (; verts.first != verts.second; ++verts.first) {
      size_t bagSize = GetProperties(*verts.first)->bag.size();
      if (bagSize > maxBag) maxBag = bagSize;
    }
    
    return maxBag - 1;
  }
  
  NTDecomp::NTDecomp() : _NTDGraph() { }
  
  NTDecomp::NTDecomp(TDecomp_p G) : _NTDGraph() {
    SetInput(G);
  }
  
  NTDecomp::NTDecomp(TDecomp_p G, TDVertex v) : _NTDGraph() {
    SetInput(G, v);
  }
  
  void NTDecomp::SetInput(TDecomp_p G) {
    TDVertex v = NULL;
    TDVertIterPair verts = G->GetVertices();
    for (; verts.first != verts.second; ++verts.first) {
      if (G->Degree(*verts.first) == 1) {
        v = *verts.first;
        //break;
      }
    }
    SetInput(G, v);
  }
  
  void NTDecomp::SetInput(TDecomp_p g, TDVertex root) {
    source_ = g;
    TDecomp_p td = g;
    PermutableGraph_p pg = td->GetSourceGraph();
    MolecularGraph_p mg = pg->GetSourceGraph();
    
    std::vector<NTDVertex> existing_vertices;
    
    // Add all vertices in TD to NTD.
    std::map<TDVertex, NTDVertex> td2ntd;
    TDVertIterPair verts = td->GetVertices();
    for (; verts.first != verts.second; ++verts.first) {
      NTDVertProp p;
      for (PermVertex v : td->GetProperties(*verts.first)->bag) {
        p.bag.insert(pg->GetProperties(v)->source);
      }
      td2ntd.emplace(*verts.first, AddVertex(p));
    }
    
    // Convert the TD to directed begining at root
    // Includes adding empty bag to every leaf
    TDVertex r = AddVertex();
    AddEdge(r, td2ntd.at(root));
    std::set<TDVertex> used;
    std::vector<TDVertex> queue;
    queue.push_back(root);
    while (queue.size()) {
      TDVertex v = queue.back();
      queue.pop_back();
      TDNbrsIterPair nbrs = td->GetNeighbours(v);
      for (; nbrs.first != nbrs.second; ++nbrs.first) {
        if (used.find(*nbrs.first) != used.end()) continue;
        AddEdge(td2ntd.at(v), td2ntd.at(*nbrs.first));
        queue.push_back(*nbrs.first);
      }
      used.insert(v);
      if (td->Degree(v) == 1 && v != root) AddEdge(td2ntd.at(v), AddVertex());
    }
    
    // Explictly add edges into bag, if not part of the TD
    bool addEdges = !Options::AssignElectrons::FPT::ADD_EDGES_TO_TD;
    if (addEdges) {
      existing_vertices.clear();
      existing_vertices.assign(GetVertices().first, GetVertices().second);
      
      for (NTDVertex leaf : existing_vertices) {
        if (Degree(leaf)) continue;
        static std::vector<MolVertPair> in_edges;
        in_edges.clear();
        
        while (InDegree(leaf)) {
          NTDVertProp* leaf_prop = GetProperties(leaf);
          static std::vector<MolVertPair> current_edges;
          current_edges.clear();
          
          for (MolVertPair uv : in_edges) {
            MolVertPair uu = std::make_pair(uv.first, uv.first);
            MolVertPair vv = std::make_pair(uv.second, uv.second);
            if (leaf_prop->bag.find(uu) == leaf_prop->bag.end()
                && leaf_prop->bag.find(vv) == leaf_prop->bag.end()) continue;
            current_edges.push_back(uv);
          }
          current_edges.swap(in_edges);
          
          for (MolVertPair uu : leaf_prop->bag) {
            if (uu.first != uu.second) continue;
            MolNbrsIterPair nbs = mg->GetNeighbours(uu.first);
            for (MolNeighboursIter n = nbs.first; n != nbs.second; ++n) {
              MolVertPair uv;
              if ((*n) > uu.first) uv = std::make_pair(uu.first, *n);
              else uv = std::make_pair(*n, uu.first);
              
              auto found = std::find(in_edges.begin(), in_edges.end(), uv);
              if (found == in_edges.end()) in_edges.push_back(uv);
            }
          }
          
          leaf_prop->bag.insert(in_edges.begin(), in_edges.end());
          leaf = *(GetPredecessors(leaf).first);
        }
      }
    }
    
    // add an intersection node between each parent/child where the child is
    // not a subset of the parent and the parent is not a subset of the child
    existing_vertices.clear();
    existing_vertices.assign(GetVertices().first, GetVertices().second);
    for (NTDVertex parent : existing_vertices) {
      NTDVertProp* p = GetProperties(parent);
      std::vector<NTDVertex> nbrs(GetNeighbours(parent).first,
                                  GetNeighbours(parent).second);
      for (NTDVertex child : nbrs) {
        NTDVertProp* c = GetProperties(child);
        std::set<MolVertPair> intersection;
        std::set_intersection(p->bag.begin(), p->bag.end(), c->bag.begin(), c->bag.end(), std::inserter(intersection, intersection.begin()));
        if (intersection == p->bag || intersection == c->bag) continue;
        NTDVertProp n;
        n.id = vert_count_++;
        n.bag.swap(intersection);
        InsertVertex(parent, child, n);
      }
    }
    
    // reduce all vertices to only have 2 children
    existing_vertices.clear();
    existing_vertices.assign(GetVertices().first, GetVertices().second);
    for (NTDVertex parent : existing_vertices) {
      while (Degree(parent) > 2) {
        NTDVertProp p;
        p.id = vert_count_++;
        p.bag = GetProperties(parent)->bag;
        
        NTDVertex new_child = AddVertex(p);
        std::vector<NTDVertex> nbrs(GetNeighbours(parent).first,
                                    GetNeighbours(parent).second);
        nbrs.pop_back();   // Must leave one existing neighbour attached
        for (NTDVertex old_child : nbrs) {
          RemoveEdge(parent, old_child);
          AddEdge(new_child, old_child);
        }
        AddEdge(parent, new_child);
        parent = new_child;
      }
    }
    
    // Reparent bad joins
    existing_vertices.clear();
    existing_vertices.assign(GetVertices().first, GetVertices().second);
    std::set<NTDVertex> removed;
    for (NTDVertex v : existing_vertices) {
      if (removed.find(v) != removed.end()) continue;
      if (Degree(v) != 2) continue;
      NTDNbrsIterPair nbrs = GetNeighbours(v);
      NTDVertex child1 = *nbrs.first;
      ++nbrs.first;
      NTDVertex child2 = *nbrs.first;
      if (Degree(child1) != 1) continue;
      if (Degree(child2) != 1) continue;
      NTDVertProp* p1 = GetProperties(child1);
      NTDVertProp* p2 = GetProperties(child2);
      if (p1->bag != p2->bag) continue;
      NTDVertex grandchild2 = *GetNeighbours(child2).first;
      AddEdge(child1, grandchild2);
      removed.emplace(child2);
      RemoveVertex(child2);
    }
    
    // add same bags after join nodes
    existing_vertices.clear();
    existing_vertices.assign(GetVertices().first, GetVertices().second);
    for (NTDVertex parent : existing_vertices) {
      if (Degree(parent) != 2) continue;
      NTDNbrsIterPair nbrs = GetNeighbours(parent);
      NTDVertex u = *nbrs.first;
      NTDVertex v = *(++nbrs.first);
      
      NTDVertProp* p_p = GetProperties(parent);
      NTDVertProp* p_u = GetProperties(u);
      NTDVertProp* p_v = GetProperties(v);
      
      if (p_u->bag != p_p->bag) {
        NTDVertProp p_n;
        p_n.id = vert_count_++;
        p_n.bag = p_p->bag;
        InsertVertex(parent, u, p_n);
      }
      
      if (p_v->bag != p_p->bag) {
        NTDVertProp p_n;
        p_n.id = vert_count_++;
        p_n.bag = p_p->bag;
        InsertVertex(parent, v, p_n);
      }
    }
    
    
    // add introduce and forget nodes as needed
    TopologicalSort(existing_vertices);
    for (NTDVertex parent : existing_vertices) {
      NTDVertProp* p = GetProperties(parent);
      std::vector<NTDVertex> nbrs(GetNeighbours(parent).first,
                                  GetNeighbours(parent).second);
      for (NTDVertex child : nbrs) {
        NTDVertProp* c = GetProperties(child);
        if (p->bag == c->bag) continue;
        std::set<MolVertPair> intersection;
        std::vector<MolVertPair> difference;
        std::set_intersection(p->bag.begin(), p->bag.end(), c->bag.begin(), c->bag.end(), std::inserter(intersection, intersection.begin()));
        if (intersection == c->bag) {
          while (p->bag.size() > c->bag.size() + 1) {
            // Introduce nodes
            difference.clear();
            std::set_difference(p->bag.begin(), p->bag.end(), c->bag.begin(), c->bag.end(), std::inserter(difference, difference.begin()));
//            size_t occurCount = 10000000;
            MolVertPair intro = difference.back();
//            for (MolVertPair v : difference) {
//              size_t currentCount = std::count(locs.begin(), locs.end(), v);
//              if (currentCount < occurCount) {
//                occurCount = currentCount;
//                intro = v;
//              }
//            }
            NTDVertProp n;
            n.id = vert_count_++;
            n.bag = c->bag;
            n.bag.insert(intro);
            NTDVertex new_child = InsertVertex(parent, child, n).first;
            child = new_child;
            c = GetProperties(child);
          }
        } else if (intersection == p->bag) {
          // Forget nodes
          while (p->bag.size() < c->bag.size() - 1) {
            difference.clear();
            std::set_difference(c->bag.begin(), c->bag.end(), p->bag.begin(), p->bag.end(), std::inserter(difference, difference.begin()));
            
            MolVertPair forget = std::make_pair(nullptr, nullptr);
            for (MolVertPair poss_forget : difference) {
              if (poss_forget.first == poss_forget.second) {
                forget = poss_forget;
                break;
              }
            }
            
            if (forget.first == nullptr) forget = difference.back();
            
            NTDVertProp n;
            n.id = vert_count_++;
            n.bag = c->bag;
            n.bag.erase(forget);
            NTDVertex new_child = InsertVertex(parent, child, n).first;
            child = new_child;
            c = GetProperties(new_child);
          }
        }
      }
    }
    
    
    
    // Identify the kind of every node
    existing_vertices.clear();
    existing_vertices.assign(GetVertices().first, GetVertices().second);
    MolVertPair nullpair = std::make_pair(nullptr, nullptr);
    for (NTDVertex v : existing_vertices) {
      NTDVertProp* p = GetProperties(v);
      std::vector<MolVertPair> difference;
      if (InDegree(v) == 0) {
        NTDVertProp* c = GetProperties(*(GetNeighbours(v).first));
        std::set_difference(c->bag.begin(), c->bag.end(), p->bag.begin(), p->bag.end(), std::inserter(difference, difference.begin()));
        p->kind = std::make_pair('R', difference.back());
      }
      else if (Degree(v) == 0) p->kind = std::make_pair('L', nullpair);
      else if (Degree(v) == 2) p->kind = std::make_pair('J', nullpair);
      else {
        NTDVertProp* c = GetProperties(*(GetNeighbours(v).first));
        if (p->bag.size() > c->bag.size()) {
          std::set_difference(p->bag.begin(), p->bag.end(), c->bag.begin(), c->bag.end(), std::inserter(difference, difference.begin()));
          p->kind = std::make_pair('I', difference.back());
        }
        else {
          std::set_difference(c->bag.begin(), c->bag.end(), p->bag.begin(), p->bag.end(), std::inserter(difference, difference.begin()));
          p->kind = std::make_pair('F', difference.back());
        }
      }
    }
    
    
    
  }
  
  /*==================*
   * Graph algorithms *
   *==================*/
  
  /*
   */
//  void NTDecomp::TopologicalSort(std::vector<NTDVertex>& sorted_out) const
//  {
//    sorted_out.clear();
//    sorted_out.reserve(NumVertices());
//    std::map<NTDVertex, int> idxMap;
//    boost::associative_property_map<std::map<NTDVertex, int>> indexMap(idxMap);
//    NTDVertIterPair verts = GetVertices();
//    for (int i = 0; verts.first != verts.second; ++verts.first, ++i) {
//      boost::put(indexMap, (*verts.first), i);
//    }
//    
//    boost::topological_sort(*graph_, std::back_inserter(sorted_out),
//                            boost::vertex_index_map(indexMap));
//  }
  
  void NTDecomp::TopologicalSort(std::vector<NTDVertex> &sorted_out) const {
    sorted_out.clear();
    sorted_out.reserve(NumVertices());
    std::vector<NTDVertex> queue;
    NTDVertIterPair vs = GetVertices();
    for (; vs.first != vs.second; ++vs.first) {
      if (InDegree(*vs.first) == 0) {
        queue.push_back(*vs.first);
        break;
      }
    }
    while (queue.size()) {
      NTDVertex v = queue.back();
      queue.pop_back();
      sorted_out.push_back(v);
      NTDNbrsIterPair nbrs = GetNeighbours(v);
      for (; nbrs.first != nbrs.second; ++nbrs.first) {
        queue.push_back(*nbrs.first);
      }
    }
    std::reverse(sorted_out.begin(), sorted_out.end());
  }
  
}
