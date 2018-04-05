/** @file counter.hpp
 *  @brief Declaration of the CountableObject template.
 *  @author Ivan Welsh
 *  @date 11 November 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_UTILS_COUNTER_HPP
#define INDIGOX_UTILS_COUNTER_HPP

#include <cstdint>

#include "../api.hpp"

namespace indigox {
  namespace utils {
    
    /** @class CountableObject counter.hpp utils/counter.hpp
     *  @brief Template class for counting instances.
     *  @tparam T class requiring counting.
     *  @details Provides a means to count instances of an object. Classes
     *  should derive from this template and call the constructor in class
     *  initalisation list.
     *  @since 0.1
     */
    template <class T>
    class CountableObject {
      
    public:
      /** @brief Default constructor.
       *  @details Sets the unique id associated with the class instance.
       */
      CountableObject() : id_(++__id_generator_) { }
      
      /** @brief Obtain the class instance unique id.
       *  @returns the unique id.
       */
      inline uid_t GetUniqueID() const { return id_; }
      
    private:
      const uid_t id_;
      static uid_t __id_generator_;
    };
    
    template <class T>
    uid_t CountableObject<T>::__id_generator_ = -1;
    
  }
}

#endif /* INDIGOX_UTILS_COUNTER_HPP */
