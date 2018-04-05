/** @file periodictable.cpp
 *  @brief Element and PeriodicTable implementation
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 7 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <cstdint>
#include <iostream>
#include <memory>
#include <sstream>
#include <vector>

#include "api.hpp"
#include "classes/periodictable.hpp"
#include "utils/common.hpp"
#include "utils/filereader.hpp"
#include "utils/options.hpp"

namespace indigox {
  
  /*
   *  Element implementation
   */
  Element::Element(uint8_t Z, std::string name, std::string symbol, float mass,
                   uint8_t group, uint8_t period, uint8_t valence,
                   uint8_t octet, uint8_t hyperOctet, float atomicR, float covR,
                   float vdwR, float chi)
  : name_(name), symbol_(symbol), grp_(group), period_(period),Z_(Z),
  val_(valence), oct_(octet), hyper_(hyperOctet), mass_(mass), radius_(atomicR),
  cov_(covR), vdw_(vdwR), chi_(chi) { }
  
  String Element::ToString() const {
    std::stringstream ss;
    ss << name_ << " (" << symbol_ << ")";
    return ss.str();
  }
  
  /*
   *  PeriodicTable implementation
   */
  
  PeriodicTable_wp PeriodicTable::instance_ = PeriodicTable_wp();
  
  /** @details If an instance does not exist, creates one and loads elemental
   *  information from the file given by the Options::PERIODIC_TABLE_FILE
   *  attribute.
   *  @returns a shared pointer to the global PeriodicTable instance.
   */
  PeriodicTable_p PeriodicTable::GetInstance() {
    PeriodicTable_p instance = instance_.lock();
    if (!instance) {
      instance.reset(new PeriodicTable());
      instance->GeneratePeriodicTable();
      instance_ = instance;
    }
    return instance;
  }
  
  /** @param z the atomic number of the element to get.
   *  @returns a shared pointer to the requested element. If no element with
   *  the given atomic number exists, the shared pointer is empty.
   */
  Element_p PeriodicTable::GetElement(uint8_t z) const {
    if (z_to_.find(z) != z_to_.end())
      return z_to_.at(z);
    return Element_p();
  }
  
  /** @param name the atomic symbol or name of the element to get.
   *  @returns a shared pointer to the requested element. If no element with
   *  the given name/symbol exists, the shared pointer is empty.
   */
  Element_p PeriodicTable::GetElement(std::string name) const {
    std::string u = utils::toUpperFirst(&name);
    if (name_to_.find(u) != name_to_.end()) {
      return name_to_.at(name).lock();
    }
    return Element_p();
  }
  
  /** @details Loads elemental data from the file pointed to by the
   *  Options::PERIODIC_TABLE_FILE attribute.
   */
  void PeriodicTable::GeneratePeriodicTable() {
    /// @todo Check for duplicate elements and skip them.
    /// @todo Ensure all columns are populated.
    std::vector<String> lines;
    if (Options::DATA_DIRECTORY.back() != '/') {
      Options::DATA_DIRECTORY.append("/");
    }
    String path = Options::DATA_DIRECTORY;
    String file = Options::PERIODIC_TABLE_FILE;
    //if (path.back() != '/') path.append("/");
    utils::FileReader fr(path + file);
    fr.GetAllItems(lines);
    
    String n = "", s = "";
    uint8_t g = 0, p = 0, z = 0, v = 0, o = 0, h = 0;
    Float m = 0.0f, r = 0.0f, c = 0.0f, V = 0.0f, e = 0.0f;
    int count = 0;
    
    for (String& line : lines) {
      switch (count) {
        case 0: z = std::stoi(line); break;
        case 1: n = utils::toUpperFirst(&line); break;
        case 2: s = utils::toUpperFirst(&line); break;
        case 3: m = std::stof(line); break;
        case 4: g = std::stoi(line); break;
        case 5: p = std::stoi(line); break;
        case 6: v = std::stoi(line); break;
        case 7: o = std::stoi(line); break;
        case 8: h = std::stoi(line); break;
        case 9: r = std::stof(line); break;
        case 10: c = std::stof(line); break;
        case 11: V = std::stof(line); break;
        case 12: e = std::stof(line); break;
        default:
          break;
      }
      
      if (count == 12) {
        Element_p elem = Element_p(new Element(z, n, s, m, g, p, v, o, h, r, c,
                                               V, e));
        z_to_.emplace(elem->GetAtomicNumber(), elem);
        name_to_.emplace(elem->GetSymbol(), elem);
        name_to_.emplace(elem->GetName(), elem);
        count = 0;
      } else count++;
    }
  }
  
  
  /*
   *  Operators
   */
  namespace ix = indigox;
  
  std::ostream& operator<<(std::ostream& s, ix::Element_p e) {
    return (s << e->ToString());
  }
  
  bool operator==(ix::Element_p l, uint8_t r) {
    return l->GetAtomicNumber() == r;
  }
  
  bool operator==(ix::Element_p l, std::string r) {
    return r.size() <= 2 ? (l->GetSymbol() == r)
    : (l->GetName() == ix::utils::toUpperFirst(&r));
  }
  
  bool operator==(ix::Element_p l, ix::Element_p r) {
    return l->GetAtomicNumber() == r->GetAtomicNumber();
  }
  
  bool operator!=(ix::Element_p l, uint8_t r) {
    return !(l == r);
  }
  
  bool operator!=(ix::Element_p l, std::string r) {
    return !(l == r);
  }
  
  bool operator!=(ix::Element_p l, ix::Element_p r) {
    return !(l == r);
  }
}
