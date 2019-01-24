//
//  electron_optimisation.cpp
//  indigox
//
//  Created by Welsh, Ivan on 12/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <algorithm>
#include <climits>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/algorithm/string.hpp>

#include <nlohmann/json.hpp>

#include "algorithm/electron_optimisation.hpp"
#include "classes/atom.hpp"
#include "classes/electron_graph.hpp"
#include "classes/molecular_graph.hpp"
#include "algorithm/formalbonds/astar.hpp"
#include "algorithm/formalbonds/fpt.hpp"
#include "algorithm/formalbonds/local_optimisation.hpp"
#include "classes/periodictable.hpp"
#include "utils/filereader.hpp"
#include "utils/options.hpp"

namespace indigox {
  enum class SortOrder {
    F_ONE,
    CL_ONE,
    BR_ONE,
    O_ONE,
    S_ONE,
    O_TWO,
    S_TWO,
    C_TWO_N_ONE,
    C_TWO_C_TWO,
    C_THREE_C_THREE,
    N_ONE,
    P_ONE,
    N_TWO,
    P_TWO,
    N_THREE,
    P_THREE,
    C_THREE_O_ONE,
    O_ONE_S_FOUR,
    C_THREE_S_ONE,
    O_ONE_P_FOUR,
    C_THREE_N_TWO,
    N_TWO_N_TWO,
    N_THREE_O_ONE,
    UNDEFINED
  };
}

using namespace indigox;
typedef Options::AssignElectrons opt_;

ElectronOpt::ElectronOpt()
: electronsToAdd_(0), molGraph_(std::make_shared<MolecularGraph>()),
elnGraph_(std::make_shared<ElectronGraph>()) { }

ElectronOpt::ElectronOpt(std::shared_ptr<MolecularGraph> G)
: electronsToAdd_(0), molGraph_(G),
elnGraph_(std::make_shared<ElectronGraph>()) { }

void ElectronOpt::SetMolecularGraph(std::shared_ptr<MolecularGraph> G) {
  molGraph_ = G;
  elnGraph_.reset(new ElectronGraph(*G));
}

size_t ElectronOpt::Run() {
  possibleLocations_.clear();
  String filename = Options::DATA_DIRECTORY + opt_::SCORE_FILE;
  if (!scores_.NumberScores() || opt_::SCORE_FILE != loaded_file) {
    scores_.LoadScoreFile(filename);
    loaded_file = String(opt_::SCORE_FILE);
  }
  
  switch (opt_::ALGORITHM) {
    case opt_::Algorithm::LOCAL_OPTIMISATION:
      algo_.reset(new algorithm::LocalOptimisation(this));
      break;
    case opt_::Algorithm::ASTAR:
      algo_.reset(new algorithm::AStarOptimisation(this));
      break;
    case opt_::Algorithm::FPT:
      algo_.reset(new algorithm::FPTOptimisation(this));
      break;
      
    default:
      throw std::runtime_error("Unsupported optimisation algorithm");
      break;
  }
  SetMolecularGraph(molGraph_);
  DetermineElectronsToAdd();
  DeterminePotentialElectronLocations();
//  for (MolVertPair& loc : possibleLocations_) {
//    MolVertex a = loc.first, b = loc.second;
//    std::cout << "(" << molGraph_->GetProperties(a)->atom->GetName() << ", " << molGraph_->GetProperties(b)->atom->GetName() << ")\n";
//  }
  if (opt_::ALGORITHM == opt_::Algorithm::ASTAR) SortPotentialLocations();
  algo_->PopulateMVP2EV();
  algo_->Run();
  size_t res = algo_->GetResultCount();
  finalScore_ = algo_->GetResultEnergy();
  return res;
}

bool ElectronOpt::ApplyElectronAssigment(Uint i) {
  return algo_->ApplyElectronAssignment(i);
}

void ElectronOpt::DetermineElectronsToAdd() {
  electronsToAdd_ = -molGraph_->GetTotalCharge();
  MolVertIterPair vs = molGraph_->GetVertices();
  for (MolVertexIter v = vs.first; v != vs.second; ++v)
    electronsToAdd_ += molGraph_->GetProperties(*v)->atom->GetElement()->GetValenceElectronCount();
  // All bonds should have order of at least 1
  // Preplace code
  ElnVertIterPair elnvs = elnGraph_->GetVertices();
  for (; elnvs.first != elnvs.second; ++elnvs.first)
    electronsToAdd_ -= elnGraph_->GetProperties(*elnvs.first)->pre_placed;
  //electronsToAdd_ -= 2 * molGraph_->NumEdges();
  // End Preplace code
  
  if (opt_::AUTO_USE_ELECTRON_PAIRS && electronsToAdd_ % 2)
    opt_::USE_ELECTRON_PAIRS = false;
  else if (opt_::AUTO_USE_ELECTRON_PAIRS)
    opt_::USE_ELECTRON_PAIRS = true;
  
  if (opt_::USE_ELECTRON_PAIRS && electronsToAdd_ % 2)
    throw std::runtime_error("Unable to handle odd number of electrons when using electron pairs.");
  else if (opt_::USE_ELECTRON_PAIRS)
    electronsToAdd_ /= 2;
}

void ElectronOpt::DeterminePotentialElectronLocations() {
  MolVertIterPair vs = molGraph_->GetVertices();
  for (MolVertexIter v = vs.first; v != vs.second; ++v) {
    int8_t octet;
    if (molGraph_->Degree(*v) > 2)
      octet = molGraph_->GetProperties(*v)->atom->GetElement()->GetHypervalentOctet();
    else
      octet = molGraph_->GetProperties(*v)->atom->GetElement()->GetOctet();
    
    int8_t bondedElectrons = 2 * molGraph_->Degree(*v);
    int8_t missingElectrons = octet - bondedElectrons;
    
    MolVertPair id = std::make_pair(*v, *v);
    // Preplace code
    ElnVertex ev = elnGraph_->GetVertex(id);
    missingElectrons -= elnGraph_->GetProperties(ev)->pre_placed;
    // End Preplace code
    
    while (missingElectrons > 0) {
      possibleLocations_.push_back(id);
      if (opt_::USE_ELECTRON_PAIRS)
        missingElectrons -= 2;
      else
        missingElectrons -= 1;
    }
  }
  
  MolEdgeIterPair es = molGraph_->GetEdges();
  for (MolEdgeIter e = es.first; e != es.second; ++e) {
    MolVertex u = molGraph_->GetSource(*e);
    MolVertex v = molGraph_->GetTarget(*e);
    if (u > v) {
      MolVertex tmp = u;
      u = v;
      v = tmp;
    }
    int8_t u_oct, v_oct, u_bond, v_bond, u_miss, v_miss;
    if (molGraph_->Degree(u) > 2)
      u_oct = molGraph_->GetProperties(u)->atom->GetElement()->GetHypervalentOctet();
    else
      u_oct = molGraph_->GetProperties(u)->atom->GetElement()->GetOctet();
    if (molGraph_->Degree(v) > 2)
      v_oct = molGraph_->GetProperties(v)->atom->GetElement()->GetHypervalentOctet();
    else
      v_oct = molGraph_->GetProperties(v)->atom->GetElement()->GetOctet();
    u_bond = 2 * molGraph_->Degree(u);
    v_bond = 2 * molGraph_->Degree(v);
    u_miss = u_oct - u_bond;
    v_miss = v_oct - v_bond;
    uint8_t order = 1;
    MolVertPair id = std::make_pair(u, v);
    // TODO: Check energy tables for available bond orders
    while (u_miss > 0 && v_miss > 0 && order <= opt_::MAXIMUM_BOND_ORDER) {
      possibleLocations_.push_back(id);
      if (!opt_::USE_ELECTRON_PAIRS)
        possibleLocations_.push_back(id);
      u_miss -= 2;
      v_miss -= 2;
      order++;
    }
  }
}

void ElectronOpt::SortPotentialLocations() {
  std::vector<ElnVertProp*> sortedUniques;
  sortedUniques.reserve(possibleLocations_.size());
  for (MolVertPair& vp : possibleLocations_) {
    ElnVertProp* p = elnGraph_->GetProperties(elnGraph_->GetVertex(vp));
    if (vp.first == vp.second) {
      MolVertProp* prop = molGraph_->GetProperties(vp.first);
      switch (prop->atom->GetElement()->GetAtomicNumber()) {
        case 7:  // Nitrogen
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::N_ONE; break;
            case 2: p->sort_score = SortOrder::N_TWO; break;
            case 3: p->sort_score = SortOrder::N_THREE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 8:  // Oxygen
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::O_ONE; break;
            case 2: p->sort_score = SortOrder::O_TWO; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 9:  // Fluorine
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::F_ONE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 15:  // Phosphorus
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::P_ONE; break;
            case 2: p->sort_score = SortOrder::P_TWO; break;
            case 3: p->sort_score = SortOrder::P_THREE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 16:  // Sulfur
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::S_ONE; break;
            case 2: p->sort_score = SortOrder::S_ONE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 17:  // Chlorine
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::CL_ONE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 35:
          switch (molGraph_->Degree(vp.first)) {
            case 1: p->sort_score = SortOrder::BR_ONE; break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        default: p->sort_score = SortOrder::UNDEFINED; break;
      }
    } else {
      MolVertProp* propa = molGraph_->GetProperties(vp.first);
      MolVertProp* propb = molGraph_->GetProperties(vp.second);
      switch (propa->atom->GetElement()->GetAtomicNumber()) {
        case 6:
          switch (propb->atom->GetElement()->GetAtomicNumber()) {
            case 6:
              if (molGraph_->Degree(vp.first) == 2
                  && molGraph_->Degree(vp.second) == 2)
                p->sort_score = SortOrder::C_TWO_C_TWO;
              else if (molGraph_->Degree(vp.first) == 3
                       && molGraph_->Degree(vp.second) == 3)
                p->sort_score = SortOrder::C_THREE_C_THREE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 7:
              if (molGraph_->Degree(vp.first) == 2
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::C_TWO_N_ONE;
              else if (molGraph_->Degree(vp.first) == 3
                       && molGraph_->Degree(vp.second) == 2)
                p->sort_score = SortOrder::C_THREE_N_TWO;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 8:
              if (molGraph_->Degree(vp.first) == 3
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::C_THREE_O_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 16:
              if (molGraph_->Degree(vp.first) == 3
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::C_THREE_S_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 7:
          switch (propb->atom->GetElement()->GetAtomicNumber()) {
            case 6:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 2)
                p->sort_score = SortOrder::C_TWO_N_ONE;
              else if (molGraph_->Degree(vp.first) == 2
                       && molGraph_->Degree(vp.second) == 3)
                p->sort_score = SortOrder::C_THREE_N_TWO;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 7:
              if (molGraph_->Degree(vp.first) == 2
                  && molGraph_->Degree(vp.second) == 2)
                p->sort_score = SortOrder::N_TWO_N_TWO;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 8:
              if (molGraph_->Degree(vp.first) == 3
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::N_THREE_O_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 8:
          switch (propb->atom->GetElement()->GetAtomicNumber()) {
            case 6:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 3)
                p->sort_score = SortOrder::C_THREE_O_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 7:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 3)
                p->sort_score = SortOrder::N_THREE_O_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 15:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 4)
                p->sort_score = SortOrder::O_ONE_P_FOUR;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 16:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 4)
                p->sort_score = SortOrder::O_ONE_S_FOUR;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 15:
          switch (propb->atom->GetElement()->GetAtomicNumber()) {
            case 8:
              if (molGraph_->Degree(vp.first) == 4
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::O_ONE_P_FOUR;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        case 16:
          switch (propb->atom->GetElement()->GetAtomicNumber()) {
            case 6:
              if (molGraph_->Degree(vp.first) == 1
                  && molGraph_->Degree(vp.second) == 3)
                p->sort_score = SortOrder::C_THREE_S_ONE;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            case 8:
              if (molGraph_->Degree(vp.first) == 4
                  && molGraph_->Degree(vp.second) == 1)
                p->sort_score = SortOrder::O_ONE_S_FOUR;
              else p->sort_score = SortOrder::UNDEFINED;
              break;
            default: p->sort_score = SortOrder::UNDEFINED; break;
          }
          break;
        default: p->sort_score = SortOrder::UNDEFINED; break;
      }
    }
    
    sortedUniques.push_back(p);
  }
  std::stable_sort(sortedUniques.begin(), sortedUniques.end(),
                   [](const ElnVertProp* lhs, const ElnVertProp* rhs) {
                     if (opt_::ALGORITHM == opt_::Algorithm::ASTAR)
                       return lhs->sort_score < rhs->sort_score;
                     else return lhs->sort_score > rhs->sort_score;
                   });
  
  for (unsigned int i = 0; i < sortedUniques.size(); ++i)
    possibleLocations_[i] = sortedUniques[i]->id;
}
