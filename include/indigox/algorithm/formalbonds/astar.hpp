//
//  astar.hpp
//  indigox
//
//  Created by Welsh, Ivan on 30/10/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef ASTAR_HPP
#define ASTAR_HPP

#include <map>
#include <queue>
#include <deque>
#include <vector>

#include <boost/dynamic_bitset.hpp>


#include "electron_optimisation_algorithm.hpp"

namespace indigox {
  namespace algorithm {
    
    typedef boost::dynamic_bitset<> VertMask;
    
    struct AStarQueueItem {
      ElnDist distribution;
      VertMask unchangeable, calculable, new_calculable;
      Score path_cost, heuristic_cost, parent_path_cost;
      uint16_t nbr_start_idx;
      uint8_t nbr_count;
    };
    
    typedef std::shared_ptr<AStarQueueItem> p_AStarQueueItem;
    
    class ItemCompare {
    public:
      bool operator()(const p_AStarQueueItem a, const p_AStarQueueItem b);
    };
    
    
    class PriorityQueue :
    public std::priority_queue<p_AStarQueueItem, std::vector<p_AStarQueueItem>,
    ItemCompare > {
    public:
      void reserve(size_t sz) { this->c.reserve(sz); }
      size_t max_size() const { return this->c.max_size(); }
      void clear() { this->c.clear(); }
    };
    
    class AStarOptimisation : public ElectronOptimisationAlgorithm {
    private:
      PriorityQueue queue_;
      
      std::vector<ElnVertex> pos2vert_;
      std::map<ElnVertex, size_t>  vert2pos_;
      
      std::vector<uint8_t> uniqueIDs_;
      std::vector<VertMask> requiredUnchangeables_;
      std::map<MolVertPair, ElnVertex> idsToVertex_;
      std::map<ElnVertex, ElnVertProp*> vertexProperties_;
      size_t vertMaskSize_;
      size_t maximumQueueSize_ = 0;
      size_t sizePerItem_;
      
    private:
      AStarOptimisation() = default;
      
    public:
      AStarOptimisation(ElectronOpt* parent);
      
    public:
      void Run() override;
      size_t GetMaximumQueueSize() { return maximumQueueSize_; }
      size_t GetSizePerQueueItem() { return sizePerItem_; }
      
    private:
      void PopulateUniqueIDs();
      void PopulateUnchangeables();
      void PopulateInitialDistribution(p_AStarQueueItem d);
      void PopulateNeighbourDistribution(p_AStarQueueItem parent, p_AStarQueueItem child);
      void Initalise();
      void DetermineCalculable(p_AStarQueueItem d);
      void CalculatePathEnergy(p_AStarQueueItem d);
      void CalculateHeuristicEnergy(p_AStarQueueItem d);
      void PromiscuousHeuristic(p_AStarQueueItem d);
      void AbstemiousHeuristic(p_AStarQueueItem d);
      void GenerateNeighbourDistributions(p_AStarQueueItem d, std::vector<ElnDist>* out_nbrs);
      
    };
  }
}

#endif /* ASTAR_HPP */
