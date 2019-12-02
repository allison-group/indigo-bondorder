//
//  java_interface.hpp
//  indigo-bondorder
//
//  Created by Welsh, Ivan on 22/11/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGO_BONDORDER_JAVA_INTERFACE_HPP
#define INDIGO_BONDORDER_JAVA_INTERFACE_HPP

#include <string>

#include "../api.hpp"

namespace indigox {
  namespace utils {
    
    String GetEliminationOrdering(String& dgf_graph);
    
  }
}

#endif /* INDIGO_BONDORDER_JAVA_INTERFACE_HPP */
