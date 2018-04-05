//
//  astar.cpp
//  indigox
//
//  Created by Welsh, Ivan on 30/10/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <algorithm>
#include <chrono>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <vector>

#include "algorithm/formalbonds/astar.hpp"

#include "utils/options.hpp"

std::ostream& operator<<(std::ostream& os, const indigox::algorithm::p_AStarQueueItem a) {
  os << "Distribution:        " << a->distribution << std::endl;
  os << "Unchangeables:       " << a->unchangeable << std::endl;
  os << "Calculable:          " << a->calculable << std::endl;
  os << "Path energy:         " << a->path_cost << std::endl;
  os << "Heuristic energy:    " << a->heuristic_cost << std::endl;
  os << "Nbr start and count: " << a->nbr_start_idx << ", " << a->nbr_count << std::endl;
  return os;
}

using namespace indigox;
using namespace algorithm;
typedef Options::AssignElectrons opt_;
typedef opt_::AStar opt_as;

bool ItemCompare::operator()(const p_AStarQueueItem a, const p_AStarQueueItem b)  {
  if (a->path_cost == opt_::INF || a->heuristic_cost == opt_::INF) {
    if (b->path_cost == opt_::INF || b->heuristic_cost == opt_::INF)
      return false;
    else
      return true;
  } else if (b->path_cost == opt_::INF || b->heuristic_cost == opt_::INF)
    return false;
  else if ((a->path_cost+a->heuristic_cost) == (b->path_cost+b->heuristic_cost))
    return (a->calculable.count()) > (b->calculable.count());
  else return (a->path_cost+a->heuristic_cost) > (b->path_cost+b->heuristic_cost);
}


AStarOptimisation::AStarOptimisation(ElectronOpt* parent)
: ElectronOptimisationAlgorithm(parent) {
  
}



void AStarOptimisation::Run() {
  Initalise();
  std::vector<ElnDist> nbrDistributions;
  nbrDistributions.reserve(8);
  
  Score targetScore = opt_::INF;
  targetScore = upperLimit_;
  
  Score initialTargetScore = targetScore;
  
  size_t count = 0;
  size_t len_limit;
  if (opt_as::MEGABYTE_LIMIT > 1048576) {
    len_limit = -1;
  } else {
    size_t test_count = 1024 * 1024;
    size_t _size = (parent_->possibleLocations_.size() / ElnDist::bits_per_block) + 1;
    _size += 2 * ((parent_->elnGraph_->NumVertices() / VertMask::bits_per_block) + 1);
    sizePerItem_ = (sizeof(AStarQueueItem) + sizeof(p_AStarQueueItem) + _size * ElnDist::bits_per_block / 4);
    test_count /= sizePerItem_;
    len_limit = opt_as::MEGABYTE_LIMIT * test_count;
  }
  
  queue_.reserve(len_limit);
  
  while (queue_.size()) {
    // Get next item
    p_AStarQueueItem source = queue_.top();
    queue_.pop();
    
    // Check exit conditions
    if (source->unchangeable.all()){
      if (targetScore == initialTargetScore)
        targetScore = source->path_cost;
      if (source->path_cost == targetScore) {
        minDistributions_.emplace_back(source->distribution, targetScore);
        if (opt_::MAXIMUM_RESULT_COUNT > 0
            && minDistributions_.size() >= opt_::MAXIMUM_RESULT_COUNT) break;
      }
      else if ((source->path_cost+source->heuristic_cost) > targetScore)
        break;
      continue;
    } else if (targetScore != initialTargetScore
               && (source->path_cost) > targetScore)
      break;
    
    // Generate neighbour queue items
    GenerateNeighbourDistributions(source, &nbrDistributions);
    for (ElnDist& item : nbrDistributions) {
      p_AStarQueueItem target = p_AStarQueueItem(new AStarQueueItem());
      target->distribution = ElnDist(item);
      
      PopulateNeighbourDistribution(source, target);
      ++count;
      
      if (target->path_cost < opt_::INF && target->heuristic_cost < opt_::INF && target->path_cost + target->heuristic_cost < initialTargetScore)
        queue_.push(target);
    }
    if (queue_.size() > maximumQueueSize_) maximumQueueSize_ = queue_.size();
    if (queue_.size() > len_limit) {
      break; // throw std::runtime_error("Exceeding A* queue memory limit");
    }
    
  }
  queue_.clear();
  if (targetScore != initialTargetScore) minScore_ = targetScore;
  else minScore_ = opt_::INF;
//  std::cout << "Maximum queue: " << maximumQueueSize_ << " items @ " << sizePerItem_ << " bytes." << std::endl;
}

void AStarOptimisation::PopulateUniqueIDs() {
  // Populate uniqueIDs_ with the uniqueIDs
  size_t num = parent_->elnGraph_->NumVertices();
  
  uniqueIDs_.reserve(num);
  
  for (unsigned int i = 0; i < num; ++i) {
    ElnVertex v = pos2vert_.at(i);
    ElnVertProp *p = parent_->elnGraph_->GetProperties(pos2vert_[i]);
    size_t count = std::count(parent_->possibleLocations_.begin(),
                              parent_->possibleLocations_.end(),
                              p->id);
    uniqueIDs_.push_back((uint8_t)count);
    idsToVertex_.emplace(p->id, v);
    vertexProperties_.emplace(v, p);
  }
}

void AStarOptimisation::PopulateUnchangeables() {
  
  requiredUnchangeables_.reserve(parent_->elnGraph_->NumVertices());
  
  for (unsigned int it = 0; it < parent_->elnGraph_->NumVertices(); ++it) {
    ElnVertex v = pos2vert_.at(it);
    ElnVertProp *p = vertexProperties_.at(v);
    VertMask unchange = VertMask(parent_->elnGraph_->NumVertices());
    unchange.reset();
    
    if (p->id.first == p->id.second) {
      // Atoms stuff
      ElnNbrsIterPair nbrs_it = parent_->elnGraph_->GetNeighbours(pos2vert_[it]);
      for (ElnNeighboursIter nbr = nbrs_it.first; nbr != nbrs_it.second; ++nbr) {
        unchange.set(vert2pos_.at(*nbr));
        ElnVertProp *nbr_p = vertexProperties_.at(*nbr);
        if (!opt_::USE_CHARGED_BOND_ENERGIES) continue;
        if (nbr_p->id.first == p->id.first)
          unchange.set(parent_->molGraph_->GetVertexIndex(nbr_p->id.second));
        else
          unchange.set(parent_->molGraph_->GetVertexIndex(nbr_p->id.first));
      }
    } else {
      unchange.set(parent_->molGraph_->GetVertexIndex(p->id.second));
      unchange.set(parent_->molGraph_->GetVertexIndex(p->id.first));
    }
    requiredUnchangeables_.push_back(unchange);
  }
}

void AStarOptimisation::Initalise() {
  uniqueIDs_.clear();
  requiredUnchangeables_.clear();
  idsToVertex_.clear();
  minDistributions_.clear();
  vertexProperties_.clear();
  pos2vert_.clear();
  vert2pos_.clear();
  queue_.clear();
  parent_->molGraph_->ResetIndicies();
  
  vertMaskSize_ = parent_->elnGraph_->NumVertices();
  
  previousDist_ = ElnDist(parent_->possibleLocations_.size());
  previousDist_.reset();
  pos2vert_.reserve(vertMaskSize_);
  ElnVertIterPair vs = parent_->elnGraph_->GetVertices();
  size_t count = 0;
  for (auto v = vs.first; v != vs.second; ++v) {
    pos2vert_.push_back(*v);
    vert2pos_.emplace(*v, count);
    ++count;
  }
  
  PopulateUniqueIDs();
  PopulateUnchangeables();
  p_AStarQueueItem init = p_AStarQueueItem(new AStarQueueItem());
  PopulateInitialDistribution(init);
  queue_.push(init);
  CalculateUpperLimit();
}

void AStarOptimisation::GenerateNeighbourDistributions(p_AStarQueueItem d, std::vector<ElnDist> *out_nbrs) {
  out_nbrs->clear();
  for (uint8_t i = 0; i <= d->nbr_count; ++i) {
    ElnDist nbr = ElnDist(d->distribution);
    for (uint8_t j = 0; j < i; ++j) {
      nbr.set(j + d->nbr_start_idx);
    }
    out_nbrs->push_back(nbr);
  }
}

void AStarOptimisation::PopulateInitialDistribution(p_AStarQueueItem d) {
  d->distribution = ElnDist(parent_->possibleLocations_.size());
  d->distribution.reset();
  
  d->unchangeable = VertMask(parent_->elnGraph_->NumVertices());
  d->unchangeable.set();
  for (unsigned int i = 0; i < d->unchangeable.size(); ++i) {
    if (uniqueIDs_.at(i) > 0)
      d->unchangeable.reset(i);
  }
  
  d->calculable = VertMask(parent_->elnGraph_->NumVertices());
  d->calculable.reset();
  DetermineCalculable(d);
  d->new_calculable = VertMask(d->calculable);
  
  d->parent_path_cost = 0;
  CalculatePathEnergy(d);
  if (d->path_cost != opt_::INF) CalculateHeuristicEnergy(d);
  
  d->nbr_start_idx = 0;
  if (parent_->possibleLocations_.size()) {
    MolVertPair startID = parent_->possibleLocations_.at(0);
    size_t startVert = vert2pos_.at(idsToVertex_.at(startID));
    d->nbr_count = (uint32_t)uniqueIDs_.at(startVert);
  } else d->nbr_count = 0;
}

void AStarOptimisation::PopulateNeighbourDistribution(p_AStarQueueItem parent, p_AStarQueueItem child) {
  
  child->unchangeable = VertMask(parent->unchangeable);
  child->calculable = VertMask(parent->calculable);
  MolVertPair startID;
  size_t startVert;
  
  if (child->distribution.count() == parent_->electronsToAdd_) {
    child->unchangeable.set();
    child->calculable.set();
    child->nbr_start_idx = (uint32_t)parent_->possibleLocations_.size();
  } else {
    startID = parent_->possibleLocations_.at(parent->nbr_start_idx);
    startVert = vert2pos_.at(idsToVertex_.at(startID));
    child->unchangeable.set(startVert);
    DetermineCalculable(child);
    child->nbr_start_idx = parent->nbr_start_idx + parent->nbr_count;
  }
  child->new_calculable = child->calculable - parent->calculable;
  
  child->parent_path_cost = parent->path_cost;
  CalculatePathEnergy(child);
  if (child->path_cost != opt_::INF) CalculateHeuristicEnergy(child);
  
  if (child->nbr_start_idx < parent_->possibleLocations_.size()) {
    startID = parent_->possibleLocations_.at(child->nbr_start_idx);
    startVert = vert2pos_.at(idsToVertex_.at(startID));
    child->nbr_count = (uint32_t)uniqueIDs_.at(startVert);
  } else child->nbr_count = 0;
  
}

void AStarOptimisation::DetermineCalculable(p_AStarQueueItem d) {
  if (d->unchangeable.all()
      || d->distribution.count() >= parent_->electronsToAdd_) {
    d->calculable.set();
    d->unchangeable.set();
    return;
  }
  
  for (unsigned int i = 0; i < d->calculable.size(); ++i) {
    if (d->calculable[i])
      continue;
    
    if (i < parent_->molGraph_->NumVertices()) {
      if (!(d->unchangeable[i]))
        continue;
      if ((d->unchangeable & requiredUnchangeables_.at(i)).count()
          != requiredUnchangeables_.at(i).count())
        continue;
      d->calculable.set(i);
    } else {
      if (!(d->unchangeable[i]))
        continue;
      if ((d->calculable & requiredUnchangeables_[i]).count() == 2)
        d->calculable.set(i);
    }
  }
}

void AStarOptimisation::CalculatePathEnergy(p_AStarQueueItem d) {
  if ((parent_->electronsToAdd_ < d->distribution.count())
      || ((d->distribution.size() - d->nbr_start_idx) < (parent_->electronsToAdd_ - d->distribution.count()))) {
    d->path_cost = opt_::INF;
  } else {
    
    SetElectronDistribution(d->distribution);
    DetermineFormalCharges();
    
    Score energy = d->parent_path_cost;
    for (size_t i = 0; i < d->calculable.size(); ++i) {
      if (d->new_calculable[i]) {
        Score v_energy = CalculateVertexEnergy(pos2vert_[i]);
        if (v_energy == opt_::INF) {
          energy = v_energy;
          break;
        } else {
          energy += v_energy;
        }
      }
    }
    
    d->path_cost = energy;
  }
}

void AStarOptimisation::CalculateHeuristicEnergy(p_AStarQueueItem d) {
  if (opt_as::HEURISTIC == opt_as::Heuristic::PROMISCUOUS)
    PromiscuousHeuristic(d);
  else if (opt_as::HEURISTIC == opt_as::Heuristic::ABSTEMIOUS)
    AbstemiousHeuristic(d);
  else
    d->heuristic_cost = 0;
}

void AStarOptimisation::PromiscuousHeuristic(p_AStarQueueItem d) {
  // Doesn't this just always return 0 or INFINITY?
  
  if (d->distribution.count() == parent_->electronsToAdd_)
    d->heuristic_cost = 0;
  else if (d->distribution.count() > parent_->electronsToAdd_)
    d->heuristic_cost = opt_::INF;
  else {
    Score h_cost = 0;
    for (unsigned int v = 0; v < d->calculable.size(); ++v) {
      ElnVertex vv = pos2vert_.at(v);
      if (d->calculable[v]) continue;
      ElnVertProp *p = vertexProperties_.at(vv);
      Score minene = opt_::INF;
      if (p->id.first == p->id.second) {  // Atom energies
        uint16_t atomic = p->atomic_number;
        
        for (uint8_t fc = 0; fc < 10; ++fc) {
          uint16_t mask = atomic + (fc << 8);
          if (parent_->scores_.find(mask) != parent_->scores_.end()) {
            Score e = parent_->scores_.at(mask);
            if (e < minene) minene = e;
          }
          
          mask += (1 << 15);
          if (parent_->scores_.find(mask) != parent_->scores_.end()) {
            Score e = parent_->scores_.at(mask);
            if (e < minene) minene = e;
          }
        }
      } else {  // Bond energies
        ElnNeighboursIter nbr = parent_->elnGraph_->GetNeighbours(vv).first;
        ElnVertProp *a = vertexProperties_.at(*nbr);
        ++nbr;
        ElnVertProp *b = vertexProperties_.at(*nbr);
        uint32_t mask_root = a->atomic_number + (b->atomic_number << 8);
        for (uint32_t a_charge = 0; a_charge < 3; ++a_charge) {
          for (uint32_t b_charge = 0; b_charge < 3; ++b_charge) {
            for (uint32_t bond_e = 1; bond_e < 9; ++bond_e) {
              uint32_t mask = mask_root + (a_charge << 16) + (b_charge << 18) + (bond_e << 20);
              if (parent_->scores_.find(mask) != parent_->scores_.end()) {
                Score e = parent_->scores_.at(mask);
                if (e < minene) minene = e;
              }
            }
          }
        }
      }
      
      
      if (minene != opt_::INF)
        h_cost += minene;
      else {
        h_cost = minene;
        break;
      }
      
    }
    d->heuristic_cost = h_cost;
  }
}

void AStarOptimisation::AbstemiousHeuristic(p_AStarQueueItem d) {
  
  if (d->distribution.count() == parent_->electronsToAdd_)
    d->heuristic_cost = 0;
  else if (d->distribution.count() > parent_->electronsToAdd_)
    d->heuristic_cost = opt_::INF;
  else {
    // Note where all the extra electrons could go
    uint32_t availableExtras = parent_->electronsToAdd_ - (uint32_t)d->distribution.count();
    std::map<ElnVertex, uint8_t> extraElectrons;
    for (unsigned int i = 0; i < uniqueIDs_.size(); ++i) {
      ElnVertex v = pos2vert_.at(i);
      if (!d->unchangeable[i] && availableExtras <= uniqueIDs_.at(i))
        extraElectrons.emplace(v, (uint8_t)availableExtras);
      else if (!d->unchangeable[i])
        extraElectrons.emplace(v, uniqueIDs_.at(i));
      else
        extraElectrons.emplace(v, 0);
    }
    std::vector<std::vector<int8_t> > attainableCharges;
    attainableCharges.reserve(d->calculable.size());
    
    Score h_cost = 0;
    
    if (opt_::USE_ELECTRON_PAIRS)
      availableExtras += availableExtras;
    
    for (unsigned int v = 0; v < d->calculable.size(); ++v) {
      ElnVertex vv = pos2vert_.at(v);
      ElnVertProp *p = vertexProperties_.at(vv);
      
      if (p->id.first == p->id.second && !d->calculable[v]) {
        // Determine how many electrons can be placed around me
        uint8_t possibleNbrElectrons = 0;
        uint8_t possibleMeElectrons = extraElectrons.at(vv);
        uint8_t step = 1;
        ElnNbrsIterPair nbr = parent_->elnGraph_->GetNeighbours(vv);
        for (ElnNeighboursIter it = nbr.first; it != nbr.second; ++it)
          possibleNbrElectrons += extraElectrons.at(*it);
        
        if (opt_::USE_ELECTRON_PAIRS)
          possibleMeElectrons += possibleMeElectrons;
        else
          possibleNbrElectrons /= 2;
        
        if ((!possibleMeElectrons % 2) && !possibleNbrElectrons)
          step = 2;
        
        uint8_t toPlace = possibleNbrElectrons + possibleMeElectrons;
        
        if (availableExtras < toPlace)
          toPlace = (uint8_t)availableExtras;
        
        // Determine possible formal charges I could attain
        std::vector<int8_t> attainable;
        attainable.reserve(toPlace+1);
        for (unsigned int i = 0; i <= toPlace; i += step)
          attainable.push_back(p->formal_charge - i);
        attainableCharges.push_back(attainable);
        
        // Determine minimum energy I can attain
        Score localMin = opt_::INF;
        uint16_t maskBase = p->atomic_number;
        for (int8_t fc : attainable) {
          uint32_t mask = maskBase + (abs(fc) << 8);
          if (fc < 0) mask += (1 << 15);
          if (parent_->scores_.find(mask) != parent_->scores_.end()
              && parent_->scores_.at(mask) < localMin)
            localMin = parent_->scores_.at(mask);
        }
        
        if (localMin == opt_::INF) {
          h_cost = localMin;
          break;
        } else {
          h_cost += localMin;
        }
        
      } else if (p->id.first == p->id.second) {
        std::vector<int8_t> attainable;
        attainable.push_back(p->formal_charge);
        attainableCharges.push_back(attainable);
      } else {
        if (d->calculable[v]) continue;
        Score localMin = opt_::INF;
        
        ElnNeighboursIter nbr = parent_->elnGraph_->GetNeighbours(vv).first;
        ElnVertProp* pa = vertexProperties_[*nbr];
        nbr++;
        ElnVertProp* pb = vertexProperties_[*nbr];
        uint32_t maskBase = pa->atomic_number + (pb->atomic_number << 8);
        if (opt_::USE_CHARGED_BOND_ENERGIES){
          if (pa->formal_charge < 0)
            maskBase += (2 << 16);
          else if (pa->formal_charge > 0)
            maskBase += (1 << 16);
          if (pb->formal_charge < 0)
            maskBase += (2 << 18);
          else if (pb->formal_charge > 0)
            maskBase += (1 << 18);
        } else if (std::count(parent_->possibleLocations_.begin(),
                              parent_->possibleLocations_.end(), p->id) == 0) {
          continue;
        }
        
        
        uint8_t toPlace = extraElectrons.at(vv);
        
        if (opt_::USE_ELECTRON_PAIRS)
          toPlace += toPlace;
        
        for (uint8_t es = 0; es <= toPlace; es += 2) {
          uint32_t mask = maskBase + ((p->electron_count + p->pre_placed + es) << 20);
          if (parent_->scores_.find(mask) != parent_->scores_.end()
              && parent_->scores_.at(mask) < localMin)
            localMin = parent_->scores_.at(mask);
          else if (opt_::USE_CHARGED_BOND_ENERGIES &&
                   parent_->scores_.find(mask & eneBitmask_) != parent_->scores_.end()
                   && parent_->scores_.at(mask & eneBitmask_) < localMin)
            localMin = parent_->scores_.at(mask & eneBitmask_);
        }
        
        if (localMin == opt_::INF) {
          h_cost = localMin;
          break;
        } else {
          h_cost += localMin;
        }
        
      }
    }
    
    d->heuristic_cost = h_cost;
  }
  
}

