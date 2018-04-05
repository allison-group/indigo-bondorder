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
  if (scores_.empty()) LoadScores();
  
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

void ElectronOpt::LoadScores() {
  PeriodicTable_p pt = PeriodicTable::GetInstance();
  if (Options::DATA_DIRECTORY.back() != '/') {
    Options::DATA_DIRECTORY.append("/");
  }
  String atmFile = Options::DATA_DIRECTORY + opt_::ATOM_ENERGY_FILE;
  utils::FileReader at(atmFile);
  String bndFile = Options::DATA_DIRECTORY + opt_::BOND_ENERGY_FILE;
  utils::FileReader bn(bndFile);
  std::vector<std::string> at_items, bn_items;
  at.GetAllItems(at_items);
  bn.GetAllItems(bn_items);
  
  int64_t global_min = LLONG_MAX;
  std::map<String, int64_t> la_mins;
  std::map<std::pair<String, String>, int64_t> lb_mins;
  
  // Convert atom energies to ints
#define NUM_DECIMAL_PLACE 5
  for (size_t i = 2; i < at_items.size(); i += 3) {
    std::string lead, tail;
    size_t sep_pos = at_items[i].find_first_of('.');
    if (sep_pos == std::string::npos)
      throw std::invalid_argument(at_items[i] + std::string(" is not a decimal value."));
    lead = at_items[i].substr(0, sep_pos);
    tail = at_items[i].substr(sep_pos + 1, NUM_DECIMAL_PLACE);
    while (tail.size() < NUM_DECIMAL_PLACE) tail.append("0");
    int64_t score;
    try {
      score = std::stoll(lead + tail);
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(at_items[i] + std::string(" is not a decimal value."));
    }
    if (la_mins.find(at_items[i-2]) == la_mins.end()) {
      la_mins.emplace(at_items[i-2], score);
    } else if (la_mins.at(at_items[i-2]) > score) {
      la_mins.at(at_items[i-2]) = score;
    }
    if (score < global_min) global_min = score;
    at_items[i] = lead + tail;
  }
  
  // Convert bond energies to ints
  for (size_t i = 3; i < bn_items.size(); i += 4) {
    std::string lead, tail;
    size_t sep_pos = bn_items[i].find_first_of('.');
    if (sep_pos == std::string::npos)
      throw std::invalid_argument(bn_items[i] + std::string(" is not a decimal value."));
    lead = bn_items[i].substr(0, sep_pos);
    tail = bn_items[i].substr(sep_pos + 1, NUM_DECIMAL_PLACE);
    while (tail.size() < NUM_DECIMAL_PLACE) tail.append("0");
    int64_t score;
    try {
      score = std::stoll(lead + tail);
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(bn_items[i] + std::string(" is not a decimal value."));
    }
    if (score < global_min) global_min = score;
    
    String atomA = bn_items[i-3];
    String atomB = bn_items[i-2];
    size_t pos_sign_a, pos_sign_b;
    pos_sign_a = atomA.find_first_of('+');
    pos_sign_b = atomB.find_first_of('+');
    
    if (pos_sign_a != std::string::npos) {
      atomA = atomA.substr(0, pos_sign_a);
    } else {
      pos_sign_a = atomA.find_first_of('-');
      if (pos_sign_a != std::string::npos) atomA = atomA.substr(0, pos_sign_a);
    }
    
    if (pos_sign_b != std::string::npos) {
      atomB = atomB.substr(0, pos_sign_b);
    } else {
      pos_sign_b = atomB.find_first_of('-');
      if (pos_sign_b != std::string::npos) atomB = atomB.substr(0, pos_sign_b);
    }
    
    std::pair<String, String> b1 = std::make_pair(atomA, atomB);
    std::pair<String, String> b2 = std::make_pair(atomB, atomA);
    if (lb_mins.find(b1) == lb_mins.end()) {
      lb_mins.emplace(b1, score);
    } else if (lb_mins.at(b1) > score) {
      lb_mins.at(b1) = score;
    }
    if (lb_mins.find(b2) == lb_mins.end()) {
      lb_mins.emplace(b2, score);
    } else if (lb_mins.at(b2) > score) {
      lb_mins.at(b2) = score;
    }
    bn_items[i] = lead + tail;
  }
  
  // Load atom energies
  for (size_t i = 0; i < at_items.size(); i += 3) {
    uint8_t Z = pt->GetElement(at_items[i])->GetAtomicNumber();
    if (Z == 0) throw std::invalid_argument(at_items[i] + std::string(" is not a valid atomic symbol."));
    int fc;
    try {
      fc = std::stoi(at_items[i + 1]);
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(at_items[i+1] + std::string(" is not an integer value."));
    }
    uint32_t k = Z + (std::abs(fc) << 8);
    if (fc < 0) k += (1 << 15);
    Score val = (Score)(std::stoll(at_items[i + 2]) - la_mins.at(at_items[i]));
    scores_.emplace(k, val);
  }
  
  // Load bond energies
  for (size_t i = 0; i < bn_items.size(); i += 4) {
    char a_sign, b_sign;
    size_t pos_sign_a, pos_sign_b;
    pos_sign_a = bn_items[i].find_first_of('+');
    pos_sign_b = bn_items[i + 1].find_first_of('+');
    
    if (pos_sign_a != std::string::npos) {
      a_sign = '+';
      bn_items[i] = bn_items[i].substr(0, pos_sign_a);
    } else {
      pos_sign_a = bn_items[i].find_first_of('-');
      if (pos_sign_a != std::string::npos) {
        a_sign = '-';
        bn_items[i] = bn_items[i].substr(0, pos_sign_a);
      } else a_sign = '0';
    }
    
    if (pos_sign_b != std::string::npos) {
      b_sign = '+';
      bn_items[i+1] = bn_items[i+1].substr(0, pos_sign_b);
    } else {
      pos_sign_b = bn_items[i+1].find_first_of('-');
      if (pos_sign_b != std::string::npos) {
        b_sign = '-';
        bn_items[i+1] = bn_items[i+1].substr(0, pos_sign_b);
      } else b_sign = '0';
    }
    
    uint8_t Za = pt->GetElement(bn_items[i])->GetAtomicNumber();
    uint8_t Zb = pt->GetElement(bn_items[i+1])->GetAtomicNumber();
    if (Za == 0) throw std::invalid_argument(at_items[i] + std::string(" is not a valid atomic symbol."));
    if (Zb == 0) throw std::invalid_argument(at_items[i+1] + std::string(" is not a valid atomic symbol."));
    
    int order;
    try {
      order = std::stoi(bn_items[i + 2]);
    } catch (const std::invalid_argument& e) {
      throw std::invalid_argument(bn_items[i+2] + std::string(" is not an integer value."));
    }
    if (order < 1) throw std::invalid_argument(bn_items[i+2] + std::string(" is an invalid bond order."));
    
    uint32_t k1 = Za + (Zb << 8) + ((order * 2) << 20);
    uint32_t k2 = Zb + (Za << 8) + ((order * 2) << 20);
    
    if (a_sign == '+') {
      k1 += (1 << 16);
      k2 += (1 << 18);
    } else if (a_sign == '-') {
      k1 += (2 << 16);
      k2 += (2 << 18);
    }
    
    if (b_sign == '+') {
      k2 += (1 << 16);
      k1 += (1 << 18);
    } else if (b_sign == '-') {
      k2 += (2 << 16);
      k1 += (2 << 18);
    }
    std::pair<String, String> bondType = std::make_pair(bn_items[i], bn_items[i+1]);
    Score val = (Score)(std::stoll(bn_items[i + 3]) - lb_mins.at(bondType));
//    Score val = (Score)(std::stoll(bn_items[i + 3]) - global_min);
    scores_.emplace(k1, val);
    if (k1 != k2) scores_.emplace(k2, val);
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
