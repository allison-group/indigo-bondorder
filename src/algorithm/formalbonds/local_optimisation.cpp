//
//  local_optimisation.cpp
//  indigox
//
//  Created by Welsh, Ivan on 13/09/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <map>
#include <set>
#include <algorithm>

#include <boost/dynamic_bitset/dynamic_bitset.hpp>

#include "algorithm/formalbonds/local_optimisation.hpp"
#include "classes/electron_graph.hpp"
#include "utils/options.hpp"

using namespace indigox;
using namespace algorithm;
typedef Options::AssignElectrons opt_;
typedef opt_::LocalOptimisation lo_;

LocalOptimisation::LocalOptimisation(ElectronOpt* parent)
: ElectronOptimisationAlgorithm(parent)
{
}

void LocalOptimisation::BuildLocationBitMasks() {
  locBitmasks_.clear();
  
  ElnDist bitmask = ElnDist(parent_->possibleLocations_.size());
  MolVertPair previousLoc = std::make_pair(nullptr, nullptr);
  
  for (size_t i = 0; i < bitmask.size(); ++i) {
    MolVertPair currentLoc = parent_->possibleLocations_[i];
    if (currentLoc == previousLoc)
      bitmask.set(i);
    else {
      if (i)
        locBitmasks_.emplace(previousLoc, bitmask);
      bitmask = ElnDist(parent_->possibleLocations_.size());
      bitmask.set(i);
      previousLoc = currentLoc;
    }
  }
  if (bitmask.size())
    locBitmasks_.emplace(previousLoc, bitmask);
}

void LocalOptimisation::Run() {
  using namespace std::chrono;
  high_resolution_clock::time_point startT, endT;
  BuildLocationBitMasks();
  std::map<ElnDist, Score> seenDists;
  std::vector<ElnDist> nbrDists;
  std::set<ElnDist> cMinDists, rMinDists, doneNbrs;
  
  previousDist_ = ElnDist(parent_->possibleLocations_.size());
  previousDist_.reset();
  
  ElnDist initDist = CalculateUpperLimit();
  cMinDists.insert(initDist);
  
  Score cMinScore = upperLimit_;
  if (cMinScore != opt_::INF) --cMinScore;
  Score rMinScore = cMinScore;
  seenDists.emplace(initDist, cMinScore);
  size_t infsEncounterd = 0;
  size_t cDistSize = cMinDists.size();
  startT = high_resolution_clock::now();
  do {
    std::set<ElnDist>::iterator begin, end;
    begin = cMinDists.begin();
    if (lo_::OPTIMISE_ALL_MINIMUMS) end = cMinDists.end();
    else end = (++cMinDists.begin());
    
    cMinScore = rMinScore;
    rMinDists.clear();
    cDistSize = cMinDists.size();
    
    while (begin != end) {
      endT = high_resolution_clock::now();
      auto t = duration_cast<milliseconds>(endT - startT).count();
      if (lo_::TIMEOUT_LIMIT != 0 && t > lo_::TIMEOUT_LIMIT) break;
//      std::set<ElnDist>::iterator doneNeighbour = doneNbrs.find(*begin);
//      if (doneNeighbour != doneNbrs.end()) {
//        ++begin;
//        continue;
//      }
//      else doneNbrs.emplace(*begin);
      GetNeighbourDistributions(*begin, &nbrDists);
      for (ElnDist n : nbrDists) {
        Score s;
        if (lo_::CACHE_RESULTS) {
          std::map<ElnDist, Score>::iterator beenSeen = seenDists.find(n);
          if (beenSeen != seenDists.end()) {
            s = beenSeen->second;
          }
          else {
            s = CalculateDistributionEnergy(n);
            if (s < opt_::INF || lo_::CACHE_INFINITIES) seenDists.emplace(n, s);
            else ++infsEncounterd;
          }
        } else {
          s = CalculateDistributionEnergy(n);
        }
        
        if (s < rMinScore) {
          rMinDists.clear();
          rMinDists.insert(n);
          rMinScore = s;
        } else if ( s == rMinScore) {
          rMinDists.insert(n);
        }
      }
      ++begin;
    }
    
    
    if (rMinScore < cMinScore) cMinDists.swap(rMinDists);
    else if (rMinScore == cMinScore) cMinDists.insert(rMinDists.begin(), rMinDists.end());
    endT = high_resolution_clock::now();
    auto t = duration_cast<milliseconds>(endT - startT).count();
    
    if (lo_::TIMEOUT_LIMIT != 0 && t > lo_::TIMEOUT_LIMIT) break;
    
  } while (cMinScore > rMinScore || cMinDists.size() != cDistSize);
  
  minDistributions_.clear();
  minDistributions_.reserve(cMinDists.size());
  if (cMinScore != Options::AssignElectrons::INF) {
    for (auto &dist : cMinDists) {
      minDistributions_.emplace_back(dist, cMinScore);
      if (opt_::MAXIMUM_RESULT_COUNT != 0
          && minDistributions_.size() >= opt_::MAXIMUM_RESULT_COUNT) break;
    }
  }
  minScore_ = cMinScore;
}

void LocalOptimisation::GetNeighbourDistributions(ElnDist dist, std::vector<ElnDist>* nbrs) {
  nbrs->clear();
  
  std::map<MolVertPair, ElnDist>::const_iterator srcBIt, srcEIt, tgtBIt, tgtEIt;
  for (srcBIt = locBitmasks_.cbegin(), srcEIt = locBitmasks_.cend(); srcBIt != srcEIt; ++srcBIt) {
    int srcIdx = (int)srcBIt->second.find_first();
    int srcLoc = srcIdx + (int)(dist & srcBIt->second).count() - 1;
    if (srcIdx > srcLoc)
      continue;
    
    for (tgtBIt = locBitmasks_.cbegin(), tgtEIt = locBitmasks_.cend(); tgtBIt != tgtEIt; ++tgtBIt) {
      if (tgtBIt == srcBIt)
        continue;
      int tgtIdx = (int)tgtBIt->second.find_first();
      int tgtCount = (int)tgtBIt->second.count();
      int tgtLoc = tgtIdx + (int)(dist & tgtBIt->second).count();
      if ((tgtCount + tgtIdx) == tgtLoc)
        continue;
      
      ElnDist newdist = ElnDist(dist);
      newdist.flip(srcLoc);
      newdist.flip(tgtLoc);
      nbrs->push_back(newdist);
    }
  }
  
}


