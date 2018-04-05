/** @file options.hpp
 *  @brief Declaration of all available options
 *  @author Ivan Welsh
 *  @date 7 January 2018
 *  @lastmodify 7 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <cstdint>
#include <set>

#include "../api.hpp"

#ifndef INDIGOX_UTILS_OPTIONS_HPP
#define INDIGOX_UTILS_OPTIONS_HPP

namespace indigox {
  
  /** @struct Options options.hpp utils/options.hpp
   *  @brief Static struct containing all runtime changeable options for
   *  indigoX.
   *  @since 0.1
   */
  struct Options {
    
    struct AssignElectrons {
      struct AStar {
        enum class Heuristic {
          PROMISCUOUS,
          ABSTEMIOUS
        };
        static Heuristic HEURISTIC;
        static uint64_t MEGABYTE_LIMIT;
        
        static void Reset();
      };
      
      struct FPT {
        enum class PermAlgo {
          RANDOM,
          QUICKBB,
          MINDEGREE,
          MINADDEDGES
        };
        static String LIBTW_JAR_FILE;
        static bool ADD_EDGES_TO_TD;
        static PermAlgo PERM_ALGO;
        static uint32_t MINIMUM_PROPAGATION_DEPTH;
        
        static void Reset();
      };
      
      struct LocalOptimisation {
        static bool OPTIMISE_ALL_MINIMUMS;
        static bool CACHE_RESULTS;
        static bool CACHE_INFINITIES;
        static uint32_t TIMEOUT_LIMIT;  // milliseconds
        
        static void Reset();
      };
      
      enum class Algorithm {
        LOCAL_OPTIMISATION,
        ASTAR,
        FPT
      };
      
      static Algorithm ALGORITHM;
      
      static String ATOM_ENERGY_FILE;
      static String BOND_ENERGY_FILE;
      
      static Score INF;
      static Uint MAXIMUM_BOND_ORDER;
      static bool USE_ELECTRON_PAIRS;
      static bool AUTO_USE_ELECTRON_PAIRS;
      static bool USE_CHARGED_BOND_ENERGIES;
      static uint32_t HIGHEST_MAGNITUDE_CHARGE;
      static bool ALLOW_CHARGED_CARBON;
      static std::set<String> ALLOWED_ELEMENTS;
      static uint32_t MAXIMUM_RESULT_COUNT;
      static bool PREPLACE_ELECTRONS;
      
      static void Reset();
    };
    
    
    
    /// @brief Default data directory path.
    static String DATA_DIRECTORY;
    
    /// @brief File in \link Options::DATA_DIRECTORY DATA_DIRECTORY \endlink
    /// containing element information for the periodic table.
    static String PERIODIC_TABLE_FILE;
    
    /// @brief Reset all attributes to their default values.
    static void Reset();
  };
  
}

#endif /* INDIGOX_UTILS_OPTIONS_HPP */
