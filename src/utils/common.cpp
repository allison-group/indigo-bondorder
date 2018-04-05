/** @file common.cpp
 *  @brief Utility functions implementation.
 *  @author Ivan Welsh
 *  @date 21 August 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <algorithm>
#include <random>

#include "api.hpp"
#include "utils/common.hpp"

/// Global namespace for all indigoX code.
namespace indigox {
  /// Nested namespace for some useful utilities.
  namespace utils {
    
    /** @param s the string to convert.
     *  @returns the upper case string.
     */
    String toUpper(const String* s){
      String t = *s;
      std::transform(t.begin(), t.end(), t.begin(), ::toupper);
      return t;
    }
    
    /** @param s the string to convert.
     *  @returns the lower case string.
     */
    String toLower(const String* s){
      String t = *s;
      std::transform(t.begin(), t.end(), t.begin(), ::tolower);
      return t;
    }
    
    /** @details Assumes the string is a single word.
     *  @param s the string to convert.
     *  @returns the lower case string with leading upper case.
     */
    String toUpperFirst(const String* s){
      String t = *s;
      std::transform(t.begin(), t.begin() + 1, t.begin(), ::toupper);
      std::transform(t.begin() + 1, t.end(), t.begin() + 1, ::tolower);
      return t;
    }
    
    /** @details Characters available are all upper and lower case letters.
     *  @param length the number of characters in the generated string.
     *  @returns a randomly generated string.
     */
    String randomString(size_t length) {
      static auto& chrs = "qwertyuiopasdfghjklzxcvbnmZAQXSWCDEVFRBGTNHYMJUKILOP";
      
      static std::mt19937 rg{std::random_device{}()};
      static std::uniform_int_distribution<size_t> pick(0, sizeof(chrs) - 2);
      
      String s;
      s.reserve(length + 4);  // no need to copy when adding extensions
      
      while(length--)
        s += chrs[pick(rg)];
      return s;
      
    }
    
  }  // namespace utils
}  // namespace indigox
