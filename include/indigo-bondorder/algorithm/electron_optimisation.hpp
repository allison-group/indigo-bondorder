//
//  electron_optimisation.hpp
//  indigo-bondorder
//
//  Created by Welsh, Ivan on 12/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ELECTRON_OPTIMISATION_HPP
#define ELECTRON_OPTIMISATION_HPP

#include <cstdint>
#include <string>
#include <vector>
#include <unordered_map>

#include "../api.hpp"

#include "../classes/electron_graph.hpp"
#include "../classes/molecular_graph.hpp"

namespace indigo_bondorder {
  namespace algorithm {
    class ElectronOptimisationAlgorithm;
    class LocalOptimisation;
    class AStarOptimisation;
    class FPTOptimisation;
  }
  
  class ElectronOpt {
    friend class indigo_bondorder::algorithm::ElectronOptimisationAlgorithm;
    friend class indigo_bondorder::algorithm::LocalOptimisation;
    friend class indigo_bondorder::algorithm::AStarOptimisation;
    friend class indigo_bondorder::algorithm::FPTOptimisation;
    
  public:
    uint32_t electronsToAdd_;
    std::vector<MolVertPair> possibleLocations_;
    std::shared_ptr<MolecularGraph> molGraph_;
    std::shared_ptr<ElectronGraph> elnGraph_;
    std::unordered_map<uint32_t, Score> scores_;
    std::shared_ptr<algorithm::ElectronOptimisationAlgorithm> algo_;
    Score finalScore_;
    
  public:
    ElectronOpt();
    ElectronOpt(std::shared_ptr<MolecularGraph> G);
    
  public:
    void SetMolecularGraph(std::shared_ptr<MolecularGraph> G);
    size_t Run();
    inline Score GetMinimisedEnergy() { return finalScore_; }
    bool ApplyElectronAssigment(Uint);
    
  private:
    void DetermineElectronsToAdd();
    void DeterminePotentialElectronLocations();
    void LoadScores();
    void SortPotentialLocations();
    
  };
}

#endif /* ELECTRON_OPTIMISATION_HPP */
