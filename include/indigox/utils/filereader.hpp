/** @file filereader.hpp
 *  @brief Declaration of the FileReader class.
 *  @author Ivan Welsh
 *  @date 13 December 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#ifndef INDIGOX_UTILS_FILEREADER_HPP
#define INDIGOX_UTILS_FILEREADER_HPP

#include <vector>

#include "../api.hpp"

namespace indigox {
  namespace utils {
    
    /** @class FileReader filereader.hpp
     *  @brief Class for reading simple text file from disk.
     *  @details Loads a simple text file from disk.
     *  @since 0.1
     */
    class FileReader {
      
    public:
      /// @brief Normal constructor.
      FileReader(const String&);
      
      /// @brief Reads the file from disk.
      void GetAllItems(std::vector<String>&);
      
    private:
      const String path_;
      
    private:
      FileReader() = default;  // No default constructor
      
    };
    
  }
}

#endif /* INDIGOX_UTILS_FILEREADER_HPP */
