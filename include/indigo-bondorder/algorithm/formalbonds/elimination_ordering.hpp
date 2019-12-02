//
//  elimination_ordering.hpp
//  indigo-bondorder
//
//  Created by Ivan Welsh on 14/01/18.
//  Copyright © 2018 Hermes Productions. All rights reserved.
//

#ifndef INDIGO_BONDORDER_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP
#define INDIGO_BONDORDER_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP

#include <vector>

#include "../../classes/permutablegraph.hpp"

namespace indigox {
  namespace algorithm {
    
    void RandomOrder(PermutableGraph_p, ElimOrder&);
    void QuickBBOrder(PermutableGraph_p, ElimOrder&);
    void MinDegreeOrder(PermutableGraph_p, ElimOrder&);
    void MinAddEdgesOrder(PermutableGraph_p, ElimOrder&);
    
  }
}

#endif /* INDIGO_BONDORDER_ALGORITHM_FORMALBONDS_ELIMINATION_ORDERING_HPP */
