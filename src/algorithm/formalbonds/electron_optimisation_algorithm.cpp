//
//  electron_optimisation_algorithm.cpp
//  indigox
//
//  Created by Welsh, Ivan on 13/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <iostream>
#include "classes/electron_graph.hpp"

#include "classes/atom.hpp"
#include "classes/bond.hpp"
#include "algorithm/electron_optimisation.hpp"
#include "algorithm/formalbonds/electron_optimisation_algorithm.hpp"
#include "utils/options.hpp"

using namespace indigox;
using namespace algorithm;

typedef Options::AssignElectrons opt_;

ElectronOptimisationAlgorithm::ElectronOptimisationAlgorithm(ElectronOpt* p)
: parent_(p)
{
  eneBitmask_ = 0;
  for (int i = 0; i < 32; ++i) {
    if (i < 16 || i > 19)
      eneBitmask_ += (1 << i);
  }
  
  // Initialise the ElectronGraph to its base state
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter v = vs.first; v != vs.second; ++v) {
    ElnVertProp* prop = parent_->elnGraph_->GetProperties(*v);
    MolVertPair id = prop->id;
    if (id.first == id.second)
      prop->electron_count = 0;
    else
      prop->electron_count = 2;
    mvp2ev_.emplace(id, *v);
  }
  
}

void ElectronOptimisationAlgorithm::PopulateMVP2EV() {
  mvp2ev_.clear();
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter v = vs.first; v != vs.second; ++v) {
    ElnVertProp* prop = parent_->elnGraph_->GetProperties(*v);
    MolVertPair id = prop->id;
    if (id.first == id.second)
      prop->electron_count = 0;
    else
      // Preplace code
      //prop->electron_count = 2;
      prop->electron_count = 0;
      // End Preplace code
    mvp2ev_.emplace(id, *v);
  }
}

bool ElectronOptimisationAlgorithm::ApplyElectronAssignment(Uint idx) {
  if (idx >= minDistributions_.size()) {
    std::cerr << "Index " << idx << " is out of range." << std::endl;
    return false;
  }
  
  std::pair<ElnDist, uint32_t> result = minDistributions_.at(idx);
  SetElectronDistribution(result.first);
  DetermineFormalCharges();
  
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter vi = vs.first; vi != vs.second; ++vi) {
    ElnVertex v = *vi;
    ElnVertProp* p = parent_->elnGraph_->GetProperties(v);
    if (p->id.first == p->id.second) {
      MolVertProp* pm = parent_->molGraph_->GetProperties(p->id.first);
      pm->atom->SetFormalCharge(p->formal_charge);
    } else if (p->id.first != p->id.second) {
      MolEdge edg = parent_->molGraph_->GetEdge(p->id.first, p->id.second).first;
      MolEdgeProp* prop = parent_->molGraph_->GetProperties(edg);
      switch (p->electron_count + p->pre_placed) {
        case 2:
          prop->bond->SetOrder(SINGLE_BOND);
          break;
        case 3:
          prop->bond->SetOrder(ONEANDAHALF_BOND);
          break;
        case 4:
          prop->bond->SetOrder(DOUBLE_BOND);
          break;
        case 5:
          prop->bond->SetOrder(TWOANDAHALF_BOND);
          break;
        case 6:
          prop->bond->SetOrder(TRIPLE_BOND);
          break;
        case 8:
          prop->bond->SetOrder(QUADRUPLE_BOND);
          break;
          
        default:
          prop->bond->SetOrder(UNDEFINED_BOND);
          break;
      }
    }
  }
  minScore_ = result.second;
  return true;
}

void ElectronOptimisationAlgorithm::SetElectronDistribution(ElnDist& dist) {
  // Add electrons per the distribution
  ElnDist changedBits = previousDist_ ^ dist;
  ElnDist lostBits = previousDist_ & changedBits;
  
  for (size_t i = 0; i < parent_->possibleLocations_.size(); ++i) {
    if (changedBits[i]) {
      ElnVertex v = mvp2ev_.at(parent_->possibleLocations_[i]);
      ElnVertProp* prop = parent_->elnGraph_->GetProperties(v);
      if (lostBits[i]) {
        if (opt_::USE_ELECTRON_PAIRS)
          prop->electron_count -= 2;
        else
          prop->electron_count--;
      } else {
        if (opt_::USE_ELECTRON_PAIRS)
          prop->electron_count += 2;
        else
          prop->electron_count++;
      }
    }
  }
  previousDist_ = dist;
}

void ElectronOptimisationAlgorithm::DetermineFormalCharges() {
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter v = vs.first; v != vs.second; ++v) {
    ElnVertProp* prop = parent_->elnGraph_->GetProperties(*v);
    MolVertPair& id = prop->id;
    
    if (id.first != id.second)
      continue;
    
    // Preplace code
    //int8_t fc = prop->valence - prop->electron_count;
    int8_t fc = prop->valence - prop->electron_count - prop->pre_placed;
    // End Preplace code
    
    ElnNbrsIterPair nbrs = parent_->elnGraph_->GetNeighbours(*v);
    uint8_t nbr_electronCount, nbr_atomicNumber;
    uint8_t atomicNumber = prop->atomic_number;
    for (ElnNeighboursIter n = nbrs.first; n != nbrs.second; ++n) {
      ElnVertProp* nbr_prop = parent_->elnGraph_->GetProperties(*n);
      // Preplace code
      //nbr_electronCount = nbr_prop->electron_count;
      nbr_electronCount = nbr_prop->electron_count + nbr_prop->pre_placed;
      // End Preplace code
      nbr_atomicNumber = nbr_prop->atomic_number;
      
      if (!(nbr_electronCount % 2))
        fc -= (nbr_electronCount / 2);
      else if (nbr_electronCount % 2
               && atomicNumber == nbr_atomicNumber
               && id.first == nbr_prop->id.first)
        fc -= ((nbr_electronCount + 1) / 2);
      else if (nbr_electronCount % 2
               && atomicNumber == nbr_atomicNumber)
        fc -= ((nbr_electronCount - 1) / 2);
      else if (nbr_electronCount % 2
               && prop->electronegativity > nbr_prop->electronegativity)
        fc -= ((nbr_electronCount + 1) / 2);
      else
        fc -= ((nbr_electronCount - 1) / 2);
    }
    
    prop->formal_charge = fc;
    
  }
}

Score ElectronOptimisationAlgorithm::CalculateDistributionEnergy(ElnDist dist) {
  SetElectronDistribution(dist);
  DetermineFormalCharges();
  Score energy = 0;
  
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter v = vs.first; v != vs.second; ++v) {
    Score v_energy = CalculateVertexEnergy(*v);
    if (v_energy == opt_::INF) {
      energy = v_energy;
      break;
    } else {
      energy += v_energy;
    }
  }
  
  return energy;
}


ElnDist ElectronOptimisationAlgorithm::CalculateUpperLimit() {
//  bool tmp_allow_charged_carbon = opt_::ALLOW_CHARGED_CARBON;
//  opt_::ALLOW_CHARGED_CARBON = true;
//  uint32_t tmp_high_charge = opt_::HIGHEST_MAGNITUDE_CHARGE;
//  opt_::HIGHEST_MAGNITUDE_CHARGE = 0;
  
  ElnDist initDist = ElnDist(previousDist_);
  initDist.reset();
  std::set<size_t> all_pos;
  std::set<size_t>::const_iterator it = all_pos.begin();
  for (size_t i = 0; i < initDist.size(); ++i) it = all_pos.insert(it, i);
  
  while (initDist.count() < parent_->electronsToAdd_) {
    Score cMin = opt_::INF;
    size_t mPos = *(all_pos.begin());
    for (size_t i : all_pos) {
      initDist.set(i);
      Score tmp = CalculateDistributionEnergy(initDist);
      if (tmp < cMin) {
        cMin = tmp;
        mPos = i;
      }
      initDist.reset(i);
    }
    initDist.set(mPos);
    all_pos.erase(mPos);
  }
  upperLimit_ = CalculateDistributionEnergy(initDist);
  
//  opt_::ALLOW_CHARGED_CARBON = tmp_allow_charged_carbon;
//  opt_::HIGHEST_MAGNITUDE_CHARGE = tmp_high_charge;
  
  if (upperLimit_ < opt_::INF) ++upperLimit_;
  return initDist;
}

Score ElectronOptimisationAlgorithm::CalculateVertexEnergy(ElnVertex& vert) {
  ElnVertProp* prop = parent_->elnGraph_->GetProperties(vert);
  MolVertPair id = prop->id;
  size_t idCount = std::count(parent_->possibleLocations_.begin(),
                              parent_->possibleLocations_.end(), id);
  
  if (id.first == id.second) { // Atom energies
    if (!opt_::ALLOW_CHARGED_CARBON && prop->atomic_number == 6
        && prop->formal_charge != 0) {
      return opt_::INF;
    }
    if (opt_::HIGHEST_MAGNITUDE_CHARGE > 0
        && (size_t)abs(prop->formal_charge) > opt_::HIGHEST_MAGNITUDE_CHARGE) {
      return opt_::INF;
    }
    uint8_t valence = prop->electron_count + prop->pre_placed;
    bool allZero = true;
    ElnNbrsIterPair nbrs = parent_->elnGraph_->GetNeighbours(vert);
    for (ElnNeighboursIter n = nbrs.first; n != nbrs.second; ++n) {
      ElnVertProp* p = parent_->elnGraph_->GetProperties(*n);
      if (std::count(parent_->possibleLocations_.begin(),
                     parent_->possibleLocations_.end(), p->id) != 0) allZero = false;
      valence += p->electron_count + p->pre_placed;
    }
    
    if (idCount == 0 && allZero) return 0;
    
    if ((parent_->elnGraph_->Degree(vert) > 2
         && valence > prop->target_hyper_octet)
        || valence > prop->target_octet){
      return opt_::INF;
    }
    
    
    uint16_t k = prop->atomic_number + (abs(prop->formal_charge) << 8);
    if (prop->formal_charge < 0)
      k += (1 << 15);
    auto pos = parent_->scores_.find(k);
    if (pos == parent_->scores_.end())
      return opt_::INF;
    else
      return pos->second;
    
  } else { // Bond energies
    if (!opt_::USE_CHARGED_BOND_ENERGIES && idCount == 0) return 0;
    ElnVertex u = parent_->elnGraph_->GetVertex(std::make_pair(id.first, id.first));
    ElnVertex v = parent_->elnGraph_->GetVertex(std::make_pair(id.second, id.second));
    ElnVertProp* u_prop = parent_->elnGraph_->GetProperties(u);
    ElnVertProp* v_prop = parent_->elnGraph_->GetProperties(v);
    
    // check valence state of both atoms
    uint8_t valence_u = u_prop->electron_count + u_prop->pre_placed;
    ElnNbrsIterPair nbrs_u = parent_->elnGraph_->GetNeighbours(u);
    for (; nbrs_u.first != nbrs_u.second; ++nbrs_u.first) {
      ElnVertProp* tmpP = parent_->elnGraph_->GetProperties(*nbrs_u.first);
      valence_u += tmpP->electron_count + tmpP->pre_placed;
    }
    if ((parent_->elnGraph_->Degree(u) > 2
         && valence_u > u_prop->target_hyper_octet)
        || valence_u > u_prop->target_octet) {
      return opt_::INF;
    }
    
    uint8_t valence_v = v_prop->electron_count + v_prop->pre_placed;
    ElnNbrsIterPair nbrs_v = parent_->elnGraph_->GetNeighbours(v);
    for (; nbrs_v.first != nbrs_v.second; ++nbrs_v.first) {
      ElnVertProp* tmpP = parent_->elnGraph_->GetProperties(*nbrs_v.first);
      valence_v += tmpP->electron_count + tmpP->pre_placed;
    }
    if ((parent_->elnGraph_->Degree(v) > 2
         && valence_v > v_prop->target_hyper_octet)
        || valence_v > v_prop->target_octet) {
      return opt_::INF;
    }
    
    uint32_t k = 0;
    k += u_prop->atomic_number;
    k += (v_prop->atomic_number << 8);
    if (opt_::USE_CHARGED_BOND_ENERGIES){
      if (u_prop->formal_charge < 0)
        k += (2 << 16);
      else if (u_prop->formal_charge > 0)
        k += (1 << 16);
      if (v_prop->formal_charge < 0)
        k += (2 << 18);
      else if (v_prop->formal_charge > 0)
        k += (1 << 18);
    }
    k += ((prop->electron_count + prop->pre_placed) << 20);
    
    auto pos = parent_->scores_.find(k);
    
    if (pos != parent_->scores_.end())
      return pos->second;
    
    pos = parent_->scores_.find(k & eneBitmask_);
    if (pos != parent_->scores_.end())
      return pos->second;
    else
      return opt_::INF;
  }
}

