//
//  local_optimisation.hpp
//  indigox
//
//  Created by Welsh, Ivan on 13/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef LOCAL_OPTIMISATION_HPP
#define LOCAL_OPTIMISATION_HPP

#include <map>
#include <vector>

#include "../../classes/molecular_graph.hpp"

#include "electron_optimisation_algorithm.hpp"

namespace indigox {
  namespace algorithm {
    
    class LocalOptimisation : public ElectronOptimisationAlgorithm {
    private:
      std::map<MolVertPair, ElnDist> locBitmasks_;
    
    private:
      LocalOptimisation() = default;
      
    public:
      LocalOptimisation(ElectronOpt* parent);
      
    public:
      void Run() override;
      
    private:
      void BuildLocationBitMasks();
      void GetNeighbourDistributions(ElnDist dist, std::vector<ElnDist>* nbrs);
      
    };
  }
}

#endif /* LOCAL_OPTIMISATION_HPP */
