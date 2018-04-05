//
//  electron_optimisation_algorithm.hpp
//  indigox
//
//  Created by Welsh, Ivan on 13/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ELECTRON_OPTIMISATION_ALGORITHM_HPP
#define ELECTRON_OPTIMISATION_ALGORITHM_HPP

#include <cstdint>
#include <map>
#include <vector>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "../../classes/electron_graph.hpp"
#include "../electron_optimisation.hpp"
#include "../../api.hpp"

namespace indigox {
  namespace algorithm {
    
    typedef boost::dynamic_bitset<> ElnDist;
    
    class ElectronOptimisationAlgorithm {
    protected:
      ElectronOpt* parent_;
      // Result storage
      std::vector<std::pair<ElnDist, uint32_t>> minDistributions_;
      std::map<MolVertPair, ElnVertex> mvp2ev_;
      Score minScore_;
      Score upperLimit_;
      
      uint32_t eneBitmask_;   // Used for checking bond energies w/wo charge
      ElnDist previousDist_;
      
    private:
      ElectronOptimisationAlgorithm() = default;
      
    public:
      ElectronOptimisationAlgorithm(ElectronOpt* parent);
      virtual ~ElectronOptimisationAlgorithm() {}
      
      //virtual void Initalise();
      virtual void Run() = 0;
      bool ApplyElectronAssignment(Uint);
      inline size_t GetResultCount() { return minDistributions_.size(); }
      inline Score GetResultEnergy() { return minScore_; }
      void PopulateMVP2EV();
      
    protected:
      void SetElectronDistribution(ElnDist& dist);
      void DetermineFormalCharges();
      ElnDist CalculateUpperLimit();
      Score CalculateVertexEnergy(ElnVertex& v);
      Score CalculateDistributionEnergy(ElnDist);
      
      
    };
    
  }
}

#endif /* ELECTRON_OPTIMISATION_ALGORITHM_HPP */
