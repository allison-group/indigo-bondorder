/** @file periodictable.hpp
 *  @brief Element and PeriodicTable declaration
 *  @author Ivan Welsh
 *  @date 5 January 2018
 *  @lastmodify 7 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_CLASSES_PERIODIC_TABLE_HPP
#define INDIGOX_CLASSES_PERIODIC_TABLE_HPP

#include <cstdint>
#include <iostream>
#include <map>
#include <memory>
#include <vector>

#include "../api.hpp"

namespace indigox {
  
  /** @class PeriodicTable periodictable.hpp classes/periodictable.hpp
   *  @brief Singleton class for storing and access elemental information.
   *  @details The PeriodicTable class provides the only means to access
   *  the Element class. Access to the instance should only be obtained
   *  using the GetInstance() method.
   *
   *  Information for the elements is stored a data file, set by the
   *  Options::PERIODIC_TABLE_FILE attribute. This attribute should be set
   *  before creating a PeriodicTable instance.
   *
   *  The file is formatted as white space seperated columns. Each line
   *  coresponds to a single element. Column order is:
   *    -# atomic number
   *    -# name of the element
   *    -# symbol for the element
   *    -# relative atomic mass in daltons
   *    -# IUPAC group number
   *    -# period
   *    -# number of valence electrons
   *    -# number of electrons in a full outer shell
   *    -# number of electrons in a full outer shell when hypervalency is allowed
   *    -# atomic radius in angstroms
   *    -# covalent radius in angstroms
   *    -# van der Waals radius in angstroms
   *    -# electronegativity on the Pauling scale
   *
   *  All columns must be populated.
   *  @since 0.1
   */
  class PeriodicTable {
  public:
    /// @brief Obtain the singleton instance of the PeriodicTable.
    static PeriodicTable_p GetInstance();
    
  public:
    /// @brief Get the element with the given atomic number.
    Element_p GetElement(uint8_t) const;
    
    /// @brief Get the element with the given name or symbol.
    Element_p GetElement(String) const;
    
    /// @returns the total number of elements in the PeriodicTable.
    inline size_t NumElements() const { return z_to_.size(); }
    
  private:
    std::map<uint8_t, Element_p> z_to_;  // owner of elements
    std::map<String, Element_wp> name_to_;
    static PeriodicTable_wp instance_;
    
  private:
    PeriodicTable() = default;
    
    /// @brief Generates the PeriodicTable data.
    void GeneratePeriodicTable();
    
  };
  
  /** @class Element periodictable.hpp classes/periodictable.hpp
   *  @brief Read only class for storing elemental information.
   *  @details Contains a large amount of relevant information pertaining
   *  to elements. No public constructors so can only be created by the
   *  PeriodicTable class.
   *  @since 0.1
   */
  class Element {
  public:
    /// @returns the relative atomic mass of the element in daltons.
    inline Float GetAtomicMass() const { return mass_; }
    
    /// @returns the atomic number of the element.
    inline uint8_t GetAtomicNumber() const { return Z_; }
    
    /// @returns the atomic radius of the element in angstroms.
    inline Float GetAtomicRadius() const { return radius_; }
    
    /// @returns the covalent radius of the element in angstroms.
    inline Float GetCovalentRadius() const { return cov_; }
    
    /// @returns the van der Waals radius of the element in angstroms.
    inline Float GetVanDerWaalsRadius() const { return vdw_; }
    
    /// @returns the name of the element.
    inline String GetName() const { return name_; }
    
    /// @returns the symbol of the element.
    inline String GetSymbol() const { return symbol_; }
    
    /// @returns the IUPAC group number the element is in.
    inline uint8_t GetGroup() const { return grp_; }
    
    /// @returns the period the element is in.
    inline uint8_t GetPeriod() const { return period_; }
    
    /// @returns the number of valence electrons the element contains.
    inline uint8_t GetValenceElectronCount() const { return val_; }
    
    /// @returns the number of electrons required for a full outer shell.
    inline uint8_t GetOctet() const { return oct_; }
    
    /// @returns the number of electrons allowed when allowing hypervalency.
    inline uint8_t GetHypervalentOctet() const { return hyper_; }
    
    /// @returns the electronegativity of the element on the Pauling scale.
    inline Float GetElectronegativity() const { return chi_; }

    /// @returns a simple string representation of the element.
    String ToString() const;
    
  private:
    const String name_, symbol_;
    const uint8_t grp_, period_, Z_, val_, oct_, hyper_;
    const Float mass_, radius_, cov_, vdw_, chi_;
    
  private:  // Only create Elements within PeriodicTable
    friend class PeriodicTable;
    Element() = default;
    Element(uint8_t, String, String, Float, uint8_t, uint8_t, uint8_t,
            uint8_t, uint8_t, Float, Float, Float, Float);
  };
  
  std::ostream& operator<<(std::ostream& s, indigox::Element_p e);
  
  bool operator==(Element_p l, uint8_t r);
  bool operator==(Element_p l, String r);
  bool operator==(Element_p l, Element_p r);
  
  bool operator!=(Element_p l, uint8_t r);
  bool operator!=(Element_p l, String r);
  bool operator!=(Element_p l, Element_p r);
}

#endif /* INDIGOX_CLASSES_PERIODIC_TABLE_HPP */
