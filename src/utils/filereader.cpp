/** @file filereader.cpp
 *  @brief Implementation of the FileReader class.
 *  @author Ivan Welsh
 *  @date 13 December 2017
 *  @lastmodify 6 January 2018
 *  @version 0.1
 *  @copyright The MIT License
 */

#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "api.hpp"
#include "utils/filereader.hpp"


namespace indigox {
  namespace utils {
    
    /** @details Sets the path of a file to read.
     *  @param file path to the file to read.
     */
    FileReader::FileReader(const String& file)
    : path_(file) {}
    
    /** @details Loads the given file and reads it itemwise (white space
     *  seperator). Lines begining with the comment character (\#) are
     *  ignored. If the comment character appears in a line, the remainder
     *  of the line is ignored.
     *  @param[out] out_items a vector to store the items from the file.
     */
    void FileReader::GetAllItems(std::vector<String> &out_items)
    {
      String next_item, line;
      std::ifstream infile;
      infile.open(path_);
      
      /// @todo Throw an exception when can't open file.
      if (!infile) {
        std::cerr << "Unable to open file: " << path_ << std::endl;
        std::exit(1);
      }
      
      while (std::getline(infile, line)) {
        std::istringstream itemiser(line);
        while (itemiser >> next_item) {
          if (next_item.at(0) == '#') break;
          out_items.push_back(next_item);
        }
      }
      
      infile.close();
    }
    
  }
}
