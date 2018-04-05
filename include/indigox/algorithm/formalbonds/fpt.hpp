//
//  fpt.hpp
//  indigox
//
//  Created by Welsh, Ivan on 4/12/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#ifndef INDIGOX_FORMALBONDS_FPT_ALGO_HPP
#define INDIGOX_FORMALBONDS_FPT_ALGO_HPP

#include <map>
#include <vector>

#include <boost/dynamic_bitset.hpp>
#include "../../api.hpp"
#include "../../classes/molecular_graph.hpp"
#include "../../classes/nicetreedecomp.hpp"

#include "electron_optimisation_algorithm.hpp"

namespace indigox {
  namespace algorithm {
    struct TDVertScore;
    
    typedef boost::dynamic_bitset<> VertMask;
    typedef VertMask ForgetMask;
    typedef VertMask BagMask;
    typedef std::multimap<Score, ForgetMask> MaskScores;
    typedef std::map<BagMask, MaskScores> BagScores;
    typedef std::map<uint32_t, BagScores> ProperScoreMatrix;
    
    class FPTOptimisation : public ElectronOptimisationAlgorithm {
      friend struct TDVertScore;
    private:
      FPTOptimisation() = default;
      
    public:
      FPTOptimisation(ElectronOpt* parent);
      
    public:
      void Run() override;
      
    private:
      void Initalise();
      void DetermineMinMax();
      void PopulateReferenceVectors();
      Score ScoreVertex(MolVertPair v, VertMask mask);
      ElnDist VertMaskToElnDist(const VertMask& m);
      
    private:
      NTDecomp_p td_;
      std::map<NTDVertex, std::shared_ptr<TDVertScore>> scorematrices_;
      std::vector<MolVertPair> possibleStates_;
      std::map<MolVertPair, VertMask> pairMasks_;
      VertMask placedMask_;
    };
    
    struct TDVertScore {
      FPTOptimisation* parent;
      NTDVertProp* tdProperties;
      ProperScoreMatrix score;
      uint32_t min_e, max_e;
      VertMask fMask;
      
      TDVertScore(FPTOptimisation* p, NTDVertProp* prop)
      : parent(p), tdProperties(prop) {}
      
      void PopulateScoreMatrix();
      VertMask IntroduceCountToMask(MolVertPair, size_t);
      void LeafPropagate();
      void ForgetPropagate(std::shared_ptr<TDVertScore> a);
      void IntroducePropagate(std::shared_ptr<TDVertScore> a);
      void JoinPropagate(std::shared_ptr<TDVertScore> a, std::shared_ptr<TDVertScore> b);
      
      String VertMaskToNiceString(VertMask);
      String ToString();
      String KindToString();
    };
    
  }
}

#endif /* INDIGOX_FORMALBONDS_FPT_ALGO_HPP */
