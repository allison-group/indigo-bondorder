//
//  electron_optimisation.hpp
//  indigox
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

#include <nlohmann/json_fwd.hpp>

#include "../api.hpp"

#include "../classes/electron_graph.hpp"
#include "../classes/molecular_graph.hpp"

namespace indigox {
  namespace algorithm {
    class ElectronOptimisationAlgorithm;
    class LocalOptimisation;
    class AStarOptimisation;
    class FPTOptimisation;
  }
  
  struct AtomMatchingType {
    using ElementMask = std::bitset<118>;
    using DegreeMask = std::bitset<16>;
    using FormalCharge = std::bitset<16>;
    
    ElementMask element;
    DegreeMask degree;
    FormalCharge formal_charge;
    
    AtomMatchingType() = default;
    AtomMatchingType(nlohmann::json& data);
    bool Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex& vert);
    
    static void SanityCheck(nlohmann::json& data);
  };
  
  struct AtomScoring {
    AtomMatchingType atom;
    std::vector<AtomMatchingType> neighbours;
    std::map<uint8_t, Score> scoring;
    Score nomatch_score;
    
    AtomScoring(nlohmann::json& data);
    bool Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex& vert);
    
    static void SanityCheck(nlohmann::json& data);
  };
  
  struct BondScoring {
    AtomMatchingType atom_a, atom_b;
    std::vector<AtomMatchingType> neighbours_a, neighbours_b;
    std::map<uint8_t, Score> scoring;
    
    BondScoring(nlohmann::json& data);
    bool Matches(std::shared_ptr<ElectronGraph> graph, ElnVertex& vert);
    
    static void SanityCheck(nlohmann::json& data);
  };
  
  struct Scoring {
    std::vector<AtomScoring> atom_scores;
    std::vector<BondScoring> bond_scores;
    
    std::map<uint8_t, std::vector<size_t>> potential_atom_matches;
    std::map<uint8_t, std::vector<size_t>> potential_bond_matches;
    
    Scoring() = default;
    
    void LoadScoreFile(std::string path);
    Score GetScore(std::shared_ptr<ElectronGraph> graph, ElnVertex& vertex);
    size_t NumberScores() {
      return atom_scores.size() + bond_scores.size();
    }
  };
  
  class ElectronOpt {
    friend class indigox::algorithm::ElectronOptimisationAlgorithm;
    friend class indigox::algorithm::LocalOptimisation;
    friend class indigox::algorithm::AStarOptimisation;
    friend class indigox::algorithm::FPTOptimisation;
    
  public:
    uint32_t electronsToAdd_;
    std::vector<MolVertPair> possibleLocations_;
    std::shared_ptr<MolecularGraph> molGraph_;
    std::shared_ptr<ElectronGraph> elnGraph_;
    static Scoring scores_;
    static String loaded_file;
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
    void SortPotentialLocations();
    
  };
}

#endif /* ELECTRON_OPTIMISATION_HPP */
