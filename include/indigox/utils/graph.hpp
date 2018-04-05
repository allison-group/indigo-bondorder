//
//  graph.hpp
//  indigox
//
//  Created by Welsh, Ivan on 10/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_UTILS_GRAPH_HPP
#define INDIGOX_UTILS_GRAPH_HPP

#include <cstdint>
#include <map>
#include <vector>

#include <boost/bimap.hpp>
#include <boost/graph/adjacency_list.hpp>

#include "../api.hpp"

namespace indigox {
  namespace utils {
    
    template <class PropType>
    struct __PropertyType {
      typedef PropType property_t;
    };
    
    struct NoProperty : __PropertyType<boost::no_property>{
    };
    
    struct DirectedGraph {
      typedef boost::bidirectionalS is_directed_t;  // So can use InDegree
    };
    
    struct UndirectedGraph {
      typedef boost::undirectedS is_directed_t;
    };
    
    // Abstracts away all boost code hopefully
    
    /* Template base class for graphs used in indigox.
     * Classes can be directed or undirected, but cannot have self loops
     * or parallel edges. Properties can be assigned to vertices and edges.
     *
     * TODO: Add handling for methods which only make sense in some graph cases.
     */
    template <
    class VertProp = NoProperty,
    class EdgeProp = NoProperty,
    class IsDirected = UndirectedGraph
    >
    class Graph {
    public:
      typedef boost::adjacency_list<
      boost::setS,
      boost::listS,
      typename IsDirected::is_directed_t,
      typename __PropertyType<VertProp>::property_t,
      typename __PropertyType<EdgeProp>::property_t > T;
      
      typedef typename boost::graph_traits<T>::vertex_descriptor VertType;
      typedef typename boost::graph_traits<T>::edge_descriptor EdgeType;
      
      typedef typename boost::graph_traits<T>::vertex_iterator VertIter;
      typedef typename boost::graph_traits<T>::edge_iterator EdgeIter;
      typedef typename boost::graph_traits<T>::adjacency_iterator NbrsIter;
      typedef typename T::inv_adjacency_iterator PredIter;
      
      typedef std::pair<VertType, VertType> VertTypePair;
      typedef std::pair<VertIter, VertIter> VertIterPair;
      typedef std::pair<EdgeIter, EdgeIter> EdgeIterPair;
      typedef std::pair<NbrsIter, NbrsIter> NbrsIterPair;
      typedef std::pair<PredIter, PredIter> PredIterPair;
      
      typedef std::pair<VertType, bool> VertBool;
      typedef std::pair<EdgeType, bool> EdgeBool;
      
      typedef boost::bimap<VertType, uid_t> VertIndexMap;
      typedef typename VertIndexMap::value_type IndexVert;
      typedef boost::bimap<EdgeType, uid_t> EdgeIndexMap;
      typedef typename EdgeIndexMap::value_type IndexEdge;
      
    protected:
      std::shared_ptr<T> graph_;
      uid_t next_edge_id_ = 0;
      uid_t next_vert_id_ = 0;
      VertIndexMap vert_idxmap_;
      EdgeIndexMap edge_idxmap_;
      
    public:
      /* Default base constructor which creates the underlaying boost graph
       * instance. Should be called in initaliser list of all classes inheriting
       * from this.
       */
      Graph() : graph_(std::make_shared<T>())
      {}
      
      Graph(const Graph& G) : graph_(std::make_shared<T>(*G.graph_))
      {}
      
      /*============================*
       * Graph modification methods *
       *============================*/
      
      /* Adds a new vertex to the graph.
       * Added vertex has no initial properties.
       * Returns the added vertex.
       */
      VertType AddVertex()
      {
        VertProp p;
        return AddVertex(p);
      }
      
      /* Adds a new vertex to the graph.
       * Added vertex has initial properties given by p.
       * Returns the added vertex.
       */
      VertType AddVertex(VertProp& p)
      {
        VertType v = boost::add_vertex(p, *graph_);
        vert_idxmap_.insert({v, next_vert_id_++});
        return v;
      }
      
      /* Adds an edge between two vertices.
       * Vertex u is used as the source and vertex v as the target.
       * Added edge has no initial properties.
       * Returns a pair of the added edge and a boolean.
       * If the boolean is false, the edge was not added and is undefined.
       */
      EdgeBool AddEdge(VertType u, VertType v)
      {
        EdgeProp p;
        return AddEdge(u, v, p);
      }
      
      /* Adds an edge between two vertices.
       * Vertex u is used as the source and vertex v as the target.
       * Added edge has initial properties given by p.
       * Returns a pair of the added edge and a boolean.
       * If the boolean is false, the edge was not added and is undefined.
       */
      EdgeBool AddEdge(VertType u, VertType v, EdgeProp p)
      {
        EdgeBool e = boost::edge(u, v, *graph_);
        if (u == v)
          return std::make_pair(e.first, false);
        if (e.second)
          return std::make_pair(e.first, false);
        e = boost::add_edge(u, v, p, *graph_);
        edge_idxmap_.insert({e.first, next_edge_id_++});
        return e;
      }
      
      /* Inserts a vertex between two existing vertices.
       * Vertices must have an existing edge between them.
       * Inserted vertex has no initial properties.
       * Returns a pair of the added vertex and a boolean.
       * If the boolean is false, the vertex was not added and is undefined.
       */
      VertBool InsertVertex(VertType u, VertType v)
      {
        VertProp p;
        return InsertVertex(u, v, p);
      }
      
      /* Inserts a vertex between two existing vertices.
       * Vertices must have an existing edge between them.
       * Existing edge is removed and new edges are added between the new vertex
       * and the existing vertices. Added edges are from u to the new vertex
       * and from the new vertex to v.
       * Inserted vertex has initial properties as given by p.
       * Added edges have properties of the existing edge.
       * Returns a pair of the added vertex and a boolean.
       * If the boolean is false, the vertex was not added and is undefined.
       */
      VertBool InsertVertex(VertType u, VertType v, VertProp p)
      {
        EdgeBool e = boost::edge(u, v, *graph_);
        if (u == v)
          return std::make_pair(u, false);
        if (!e.second)
          return std::make_pair(v, false);
        EdgeProp ep = (*graph_)[e.first];
        VertType i = AddVertex(p);
        RemoveEdge(e.first);
        if (AddEdge(u, i, ep).second && AddEdge(i, v, ep).second)
        {
          return std::make_pair(i, true);
        }
        RemoveVertex(i);
        return std::make_pair(i, false);
      }
      
      /* Removes a vertex.
       * Vertex v is removed from the graph.
       * If v does not exist in the graph, removal silently fails.
       */
      void RemoveVertex(VertType v)
      {
        boost::clear_vertex(v, *graph_);
        boost::remove_vertex(v, *graph_);
        if (vert_idxmap_.left.find(v) != vert_idxmap_.left.end()) {
          vert_idxmap_.left.erase(v);
        }
        
      }
      
      /* Removes an edge.
       * Edge e is removed from the graph.
       * If e does not exist in the graph, removal silently fails.
       */
      void RemoveEdge(EdgeType e)
      {
        VertType u = boost::source(e, *graph_);
        VertType v = boost::target(e, *graph_);
        RemoveEdge(u, v);
      }
      
      /* Removes an edge between two vertices.
       * If an edge exists between u and v, it is removed from the graph.
       * If an edge does not exist between u and v, removal silently fails.
       */
      void RemoveEdge(VertType u, VertType v)
      {
        EdgeBool e = GetEdge(u, v);
        if (!e.second) return;
        if (edge_idxmap_.left.find(e.first) != edge_idxmap_.left.end()) {
          edge_idxmap_.left.erase(e.first);
        }
        boost::remove_edge(u, v, *graph_);
      }
      
      /* Removes all edges and vertices from the graph
       */
      void Clear()
      {
        graph_->clear();
        vert_idxmap_.clear();
        edge_idxmap_.clear();
        next_vert_id_ = 0;
        next_edge_id_ = 0;
      }
      
      /*===============================*
       * Graph access and data methods *
       *===============================*/
      
      /* Return the number of vertices in the graph.
       */
      size_t NumVertices() const
      {
        return boost::num_vertices(*graph_);
      }
      
      /* Return the number of edges in the graph.
       */
      size_t NumEdges() const
      {
        return boost::num_edges(*graph_);
      }
      
      /* Return the degree of a vertex.
       * The degree of a vertex is the number of edges from it.
       */
      size_t Degree(VertType v) const
      {
        return boost::out_degree(v, *graph_);
      }
      
      
      size_t InDegree(VertType v) const
      {
        return boost::in_degree(v, *graph_);
      }
      
      /* Return a pair of iterators over the neighbours of a vertex.
       * The first pair member is the start of the iteration range, the second is the end.
       */
      NbrsIterPair GetNeighbours(VertType v) const
      {
        return boost::adjacent_vertices(v, *graph_);
      }
      
      /* Return a pair of iterators over the neighbours incident to a vertex.
       * The first pair member is the start of the iteration range, the second is the end.
       */
      PredIterPair GetPredecessors(VertType v) const
      {
        return boost::inv_adjacent_vertices(v, *graph_);
      }
      
      /* Return the pair of vertices which an edge is between.
       * The first pair member is the source, second is the
       */
      VertTypePair GetVertices(EdgeType e) const
      {
        VertType u = boost::source(e, *graph_);
        VertType v = boost::target(e, *graph_);
        return std::make_pair(u, v);
      }
      
      /* Return a pair of iterators over the vertices in the graph.
       * First pair member is the start of the iteration range.
       * Second pair member is the end of the iteration range.
       */
      VertIterPair GetVertices() const
      {
        return boost::vertices(*graph_);
      }
      
      /* Return a pair of iterators over the edges in the graph.
       * First pair member is the start of the iteration range.
       * Second pair member is the end of the iteration range.
       */
      EdgeIterPair GetEdges() const
      {
        return boost::edges(*graph_);
      }
      
      /* Return a pointer to the properties of vertex v.
       * No checking is performed to ensure the vertex exists in the graph.
       */
      VertProp* GetProperties(VertType v) const
      {
        return &(*graph_)[v];
      }
      
      /* Return a pointer to the properties of edge e.
       * No checking is performed to ensure the edge exists in the graph.
       */
      EdgeProp* GetProperties(EdgeType e) const
      {
        return &(*graph_)[e];
      }
      
      /* Return the edge between two vertices.
       * Edge is returned as an edge-bool pair.
       * First pair member is the edge, second is the boolean.
       * If the boolean is false, there is no edge between the
       * two vertices and it is undefined.
       */
      EdgeBool GetEdge(VertType u, VertType v) const
      {
        return boost::edge(u, v, *graph_);
      }
      
      EdgeBool GetEdgeByIndex(uid_t idx) const
      {
        if (edge_idxmap_.right.find(idx) != edge_idxmap_.right.end()) {
          return EdgeBool(edge_idxmap_.right.at(idx), true);
        }
        return EdgeBool(EdgeType(), false);
      }
      
      /* Return the source vertex of an edge.
       */
      VertType GetSource(EdgeType e) const
      {
        return boost::source(e, *graph_);
      }
      
      /* Return the target vertex of an edge.
       */
      VertType GetTarget(EdgeType e) const
      {
        return boost::target(e, *graph_);
      }
      
      VertBool GetVertexByIndex(uid_t idx) const {
        if (vert_idxmap_.right.find(idx) != vert_idxmap_.right.end()) {
          return VertBool(vert_idxmap_.right.at(idx), true);
        }
        return VertBool(VertType(), false);
      }
      
      uid_t GetVertexIndex(VertType v) const {
        if (vert_idxmap_.left.find(v) != vert_idxmap_.left.end()) {
          return vert_idxmap_.left.at(v);
        }
        return -1;
      }
      
      uid_t GetEdgeIndex(EdgeType e) const {
        if (edge_idxmap_.left.find(e) != edge_idxmap_.left.end()) {
          return edge_idxmap_.left.at(e);
        }
        return -1;
      }
      
      void ResetIndicies() {
        vert_idxmap_.clear();
        edge_idxmap_.clear();
        next_vert_id_ = 0;
        next_edge_id_ = 0;
        
        VertIterPair vs = GetVertices();
        while (vs.first != vs.second) {
          vert_idxmap_.insert({*vs.first, next_vert_id_++});
          ++vs.first;
        }
        
        EdgeIterPair es = GetEdges();
        while (es.first != es.second) {
          edge_idxmap_.insert({*es.first, next_edge_id_++});
          ++es.first;
        }
      }
      
    };
  }
}

#endif /* INDIGOX_UTILS_GRAPH_HPP */

