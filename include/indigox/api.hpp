//
//  api.hpp
//  indigox
//
//  Created by Ivan Welsh on 7/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//
#include <cstdint>
#include <memory>
#include <string>

#ifndef INDIGOX_API_HPP
#define INDIGOX_API_HPP

namespace indigox {
  // General typedefs
  typedef int Int;
  typedef unsigned int Uint;
  typedef float Float;
  typedef std::string String;
  typedef unsigned long uid_t;
  typedef unsigned int Score;
  
  // All simple classes
  class Atom;
  class Bond;
  class MolecularGraph;
  class Molecule;
  class Element;
  class NTDecomp;
  class PeriodicTable;
  class PermutableGraph;
  class TDecomp;
  
  // shared_ptr of simple classes
  typedef std::shared_ptr<Atom> Atom_p;
  typedef std::shared_ptr<Bond> Bond_p;
  typedef std::shared_ptr<MolecularGraph> MolecularGraph_p;
  typedef std::shared_ptr<Molecule> Molecule_p;
  typedef std::shared_ptr<Element> Element_p;
  typedef std::shared_ptr<PeriodicTable> PeriodicTable_p;
  typedef std::shared_ptr<PermutableGraph> PermutableGraph_p;
  typedef std::shared_ptr<TDecomp> TDecomp_p;
  typedef std::shared_ptr<NTDecomp> NTDecomp_p;
  
  // weak_ptr of simple classes
  typedef std::weak_ptr<Atom> Atom_wp;
  typedef std::weak_ptr<Bond> Bond_wp;
  typedef std::weak_ptr<Molecule> Molecule_wp;
  typedef std::weak_ptr<MolecularGraph> MolecularGraph_wp;
  typedef std::weak_ptr<Element> Element_wp;
  typedef std::weak_ptr<PeriodicTable> PeriodicTable_wp;
  
  namespace utils {
    // CountableObject related
    template <class T>
    class CountableObject;
    
    // Graph related
    template <class PropType>
    struct __PropertyType;
    struct NoProperty;
    struct DirectedGraph;
    struct UndirectedGraph;
    template <class VertProp, class EdgeProp, class UndirectedGraph>
    class Graph;
  }
  
}


#endif /* INDIGOX_API_HPP */

