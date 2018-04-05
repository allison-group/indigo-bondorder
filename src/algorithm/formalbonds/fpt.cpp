//
//  fpt.cpp
//  indigox
//
//  Created by Welsh, Ivan on 4/12/17.
//  Copyright © 2017 Allison Group. All rights reserved.
//

#include <algorithm>
#include <iostream>
#include <chrono>
#include <vector>

#include "algorithm/formalbonds/fpt.hpp"
#include "algorithm/formalbonds/elimination_ordering.hpp"
#include "utils/options.hpp"
#include "classes/nicetreedecomp.hpp"
#include "classes/treedecomp.hpp"
#include "classes/permutablegraph.hpp"

using namespace indigox;
using namespace algorithm;

typedef Options::AssignElectrons opt_;

template <class T>
struct InnerIterators {
  typename T::const_iterator begin;
  typename T::const_iterator end;
  typename T::const_iterator current;
};

template <class T, class ICT = std::vector<T>, class OCT = std::vector<ICT>>
class CartesianProduct {
public:
  typedef T type;
  typedef ICT innerType;
  typedef OCT outerType;
  
  typedef typename innerType::const_iterator innerIter;
  typedef typename outerType::const_iterator outerIter;
  typedef InnerIterators<innerType> innerIters;
  typedef typename std::vector<innerIters> innerItersC;
  
private:
  innerItersC iters_;
  bool atEnd_;
  
private:
  CartesianProduct() = default;
  
public:
  CartesianProduct(outerIter begin, outerIter end)
  : atEnd_(false)
  {
    for (; begin != end; ++begin) {
      innerIters it = {(*begin).begin(), (*begin).end(), (*begin).begin()};
      iters_.push_back(it);
    }
    if (!iters_.size()) atEnd_ = true;
  }
  
  CartesianProduct(innerType& a, innerType& b)
  : atEnd_(false)
  {
    innerIters ita = {a.begin(), a.end(), a.begin()};
    innerIters itb = {b.begin(), b.end(), b.begin()};
    iters_.push_back(ita);
    iters_.push_back(itb);
    if (!iters_.size()) atEnd_ = true;
    if (a.begin() == a.end() || b.begin() == b.end()) atEnd_ = true;
  }
  
  void GetNextProduct(innerType& c)
  {
    c.clear();
    if (atEnd_) return;
    
    for (auto it = iters_.begin(); it != iters_.end(); ++it)
      c.push_back(*(it->current));
    
    for (auto it = iters_.begin(); ; ) {
      ++(it->current);
      if (it->current == it->end) {
        if (it+1 == iters_.end()) atEnd_ = true;
        else {
          it->current = it->begin;
          ++it;
        }
      } else break;
    }
    return;
    
  }
  
};

FPTOptimisation::FPTOptimisation(ElectronOpt* parent)
: ElectronOptimisationAlgorithm(parent)
{ }

void FPTOptimisation::Run()
{
#ifdef DEBUG_DATA
  std::chrono::high_resolution_clock::time_point start, end;
  start = std::chrono::high_resolution_clock::now();
#endif
  Initalise();
  std::vector<NTDVertex> vertOrder;
  vertOrder.reserve(td_->NumVertices());
  td_->TopologicalSort(vertOrder);
  std::map<NTDVertex, std::vector<NTDVertex>> deleteAfters;
  
  int count = 0;
  for (NTDVertex v : vertOrder) {
    NTDVertex parent = *(td_->GetPredecessors(v).first);
    if (deleteAfters.find(parent) == deleteAfters.end()) {
      std::vector<NTDVertex> todels;
      deleteAfters.emplace(parent, todels);
    }
    count++;
    deleteAfters.at(parent).push_back(v);
  }
  
  for (NTDVertex v : vertOrder) {
    std::shared_ptr<TDVertScore> s = scorematrices_.at(v);
    NTDVertProp* p = td_->GetProperties(v);
    NTDNbrsIterPair nbrs = td_->GetNeighbours(v);
    NTDNeighboursIter n = nbrs.first;
    std::shared_ptr<TDVertScore> csa, csb;
//    std::cout << td_->GetVertexIndex(v) << " " << s->KindToString() << std::endl;
//    std::cout << "PrePropagate:" << std::endl << s->ToString();
    switch (p->kind.first) {
      case 'L':
        s->LeafPropagate();
        break;
      case 'F':
        csa = scorematrices_.at(*n);
//        std::cout << "Previous " << td_->GetVertexIndex(*n) << std::endl;
        s->ForgetPropagate(csa);
        break;
      case 'I':
        csa = scorematrices_.at(*n);
//        std::cout << "Previous " << td_->GetVertexIndex(*n) << std::endl;
        s->IntroducePropagate(csa);
        break;
      case 'J':
        csa = scorematrices_.at(*n);
//        std::cout << "Joining " << td_->GetVertexIndex(*n) << " and ";
        n++;
        csb = scorematrices_.at(*n);
//        std::cout << td_->GetVertexIndex(*n) << std::endl;
        s->JoinPropagate(csa, csb);
        break;
      case 'R':
        csa = scorematrices_.at(*n);
//        std::cout << "Previous " << td_->GetVertexIndex(*n) << std::endl;
        s->ForgetPropagate(csa);
        // Is the last node, so output the energies found?
        minScore_ = opt_::INF;
        for (auto& e : s->score) {
          for (auto& bm : e.second) {
            for (auto& sm : bm.second) {
              if (sm.first < minScore_) {
                minScore_ = sm.first;
                minDistributions_.clear();
                minDistributions_.emplace_back(VertMaskToElnDist(sm.second), sm.first);
              } else if (sm.first == minScore_) {
                minDistributions_.emplace_back(VertMaskToElnDist(sm.second), sm.first);
              }
            }
          }
        }
        break;
    }
    
    if (deleteAfters.find(v) != deleteAfters.end()) {
      for (NTDVertex u : deleteAfters.at(v))
        scorematrices_.at(u).reset();
    }
    
  }
#ifdef DEBUG_DATA
  end = std::chrono::high_resolution_clock::now();
  auto micro = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "FPT: " << std::endl;
  std::cout << "Run time: " << micro.count() << " microseconds." << std::endl;
  std::cout << "Min score: " << minScore_ << std::endl;
  std::cout << "Minimum count: " << minDistributions_.size() << std::endl;
  std::cout << std::endl;
#endif
  
}

void FPTOptimisation::Initalise()
{
  possibleStates_.clear();
  scorematrices_.clear();
  pairMasks_.clear();
  td_.reset();
  
  previousDist_ = ElnDist(parent_->possibleLocations_.size());
  previousDist_.reset();
  PermutableGraph_p PG = PermutableGraph_p(new PermutableGraph(parent_->molGraph_));
  ElimOrder order;
  switch (opt_::FPT::PERM_ALGO) {
    case opt_::FPT::PermAlgo::RANDOM:
      RandomOrder(PG, order);
      break;
    case opt_::FPT::PermAlgo::QUICKBB:
      QuickBBOrder(PG, order);
      break;
    case opt_::FPT::PermAlgo::MINDEGREE:
      MinDegreeOrder(PG, order);
      break;
    case opt_::FPT::PermAlgo::MINADDEDGES:
      MinAddEdgesOrder(PG, order);
      break;
  }
  
  TDecomp_p TD = TDecomp_p(new TDecomp(PG, order));
  td_ = NTDecomp_p(new NTDecomp(TD));
  PopulateReferenceVectors();
  DetermineMinMax();
  CalculateUpperLimit();
}

// Populates the all_forgets vector which lists all MolVertPairs which
// can be forgotten, multiple times for those that can have electrons
// added.
void FPTOptimisation::PopulateReferenceVectors()
{
  size_t targetSize = parent_->elnGraph_->NumVertices() +
  parent_->possibleLocations_.size();
  possibleStates_.reserve(targetSize);
  placedMask_ = VertMask(targetSize);
  
  ElnVertIterPair ids = parent_->elnGraph_->GetVertices();
  for (ElnVertexIter id = ids.first; id != ids.second; ++id) {
    ElnVertProp* p = parent_->elnGraph_->GetProperties(*id);
    MolVertPair vp = p->id;
    VertMask mask(targetSize);
    
    size_t count = std::count(parent_->possibleLocations_.begin(),
                              parent_->possibleLocations_.end(),
                              vp) + 1;
    for (size_t i = 0; i < count; ++i)
      mask.set(i + possibleStates_.size());
    for (size_t i = 0; count > 1 && i < count - 1; ++i)
      placedMask_.set(i + possibleStates_.size());
    
    pairMasks_.emplace(vp, mask);
    while (count > 0) {
      possibleStates_.push_back(vp);
      count--;
    }
  }
}

// Determines the minimum and maximum numbers of electrons to be forgotten
// in the TD below, and including, every vertex
void FPTOptimisation::DetermineMinMax()
{
  NTDVertIterPair vs = td_->GetVertices();
  for (NTDVertexIter v = vs.first; v != vs.second; ++v) {
    NTDVertProp* p = td_->GetProperties(*v);
    std::shared_ptr<TDVertScore> score(new TDVertScore(this, p));
    std::vector<MolVertPair> forgottens;
    std::vector<NTDVertex> descendants;
    descendants.push_back(*v);
    while (descendants.size()) {
      NTDVertex u = descendants.back();
      descendants.pop_back();
      p = td_->GetProperties(u);
      if (p->kind.first == 'F' || p->kind.first == 'R')
        forgottens.push_back(p->kind.second);
      NTDNbrsIterPair nbsp = td_->GetNeighbours(u);
      for (NTDNeighboursIter n = nbsp.first; n != nbsp.second; ++n)
        descendants.push_back(*n);
    }
    
    int n_pos = 0;
    for (MolVertPair id : forgottens)
      n_pos += std::count(parent_->possibleLocations_.begin(),
                          parent_->possibleLocations_.end(),
                          id);
    int n_pos_invert = (int)parent_->possibleLocations_.size() - n_pos;
    score->min_e = std::max(0, (int)parent_->electronsToAdd_ - n_pos_invert);
    score->max_e = std::min((int)parent_->electronsToAdd_, n_pos);
    scorematrices_.emplace(*v, score);
  }
}

VertMask TDVertScore::IntroduceCountToMask(MolVertPair v, size_t count)
{
  VertMask mask = VertMask(parent->pairMasks_.at(v));
  while (mask.count() > count) mask.reset(mask.find_first());
  return mask;
}


void TDVertScore::ForgetPropagate(std::shared_ptr<TDVertScore> a)
{
  MolVertPair f = tdProperties->kind.second;
  fMask = a->fMask | parent->pairMasks_.at(f);
  for (uint32_t e = min_e; e <= max_e; ++e) {
    BagScores tmpBag;
    score.emplace(e, tmpBag);
  }
  
  VertMask tmpF = parent->placedMask_ & fMask;
  for (auto& e : a->score) {
    for (auto& bm : e.second) {
      BagMask currentBagMask = bm.first - fMask;
//      std::cout << "Forget jobs @" << a->min_e << "-" << e.first << "-" << a->max_e << ": " << bm.second.size() << std::endl;
      for (auto& maskScore : bm.second) {
        VertMask checkMask = maskScore.second & tmpF;
        uint32_t eNew = (uint32_t)checkMask.count();
        if (eNew < min_e || eNew > max_e) continue;
        Score s = parent->ScoreVertex(f, maskScore.second);
        if (s == opt_::INF) continue;
        s += maskScore.first;
        if (s > parent->upperLimit_) continue;
        if (score.at(eNew).find(currentBagMask) == score.at(eNew).end()) {
          MaskScores tmpScores;
          score.at(eNew).emplace(currentBagMask, tmpScores);
        }
        if (opt_::MAXIMUM_RESULT_COUNT == 0
            || score.at(eNew).at(currentBagMask).count(s) < opt_::MAXIMUM_RESULT_COUNT) {
          score.at(eNew).at(currentBagMask).emplace(s, maskScore.second);
        }
      }
    }
  }
  
  if (opt_::FPT::MINIMUM_PROPAGATION_DEPTH > 0) {
    for (auto& e : score) {
      for (auto& bm : e.second) {
        std::set<Score> minScores;
        for (auto& ms : bm.second) {
          if (minScores.size() < opt_::FPT::MINIMUM_PROPAGATION_DEPTH) {
            minScores.emplace(ms.first);
          } else {
            Score maxScore = *std::max_element(minScores.begin(), minScores.end());
            if (ms.first < maxScore) minScores.emplace(ms.first);
            if (minScores.size() > opt_::FPT::MINIMUM_PROPAGATION_DEPTH)
              minScores.erase(maxScore);
          }
        }
        Score maxScore = *std::max_element(minScores.begin(), minScores.end());
        for (auto it = bm.second.begin(); it != bm.second.end();) {
          if (it->first > maxScore) it = bm.second.erase(it);
          else ++it;
        }
      }
    }
  }
  
//  if (opt_::MAXIMUM_RESULT_COUNT > 0) {
//    uint32_t x = opt_::MAXIMUM_RESULT_COUNT;
//    for (auto& e : score) {
//      for (auto& bm : e.second) {
//        std::multiset<Score> seenScores;
//        for (auto it = bm.second.begin(); it != bm.second.end();) {
//          seenScores.emplace(it->second);
//          if (seenScores.count(it->second) > x) {
//            it = bm.second.erase(it);
//          } else {
//            ++it;
//          }
//        }
//      }
//    }
//  }
  
  a->score.clear();
}

void TDVertScore::LeafPropagate() {
  fMask = VertMask(parent->placedMask_.size());
  MaskScores tmpScores;
  BagScores tmpBags;
  score.emplace(0,tmpBags);
  score.at(0).emplace(fMask, tmpScores);
  score.at(0).at(fMask).emplace(0, fMask);
}

void TDVertScore::IntroducePropagate(std::shared_ptr<TDVertScore> a)
{
  MolVertPair introduced = tdProperties->kind.second;
  fMask = a->fMask;
  size_t count = parent->pairMasks_.at(introduced).count();
  for (uint32_t e = min_e; e <= max_e; ++e) {
    BagScores tmpBags;
    score.emplace(e, tmpBags);
  }
  
  for (size_t i = 1; i <= count; ++i) {
    BagMask iMask = IntroduceCountToMask(introduced, i);
    for (uint32_t e = min_e; e <= max_e; ++e) {
      for (auto& bmask : a->score.at(e)) {
        BagMask newMask = bmask.first | iMask;
        MaskScores tmpScores;
        score.at(e).emplace(newMask, tmpScores);
        for (auto& ms : bmask.second) {
          score.at(e).at(newMask).emplace(ms.first, ms.second | iMask);
        }
      }
    }
  }
  a->score.clear();
}

void TDVertScore::JoinPropagate(std::shared_ptr<TDVertScore> a, std::shared_ptr<TDVertScore> b)
{
  fMask = a->fMask | b->fMask;
  VertMask notFMask = VertMask(fMask);
  notFMask.flip();
  for (uint32_t e = min_e; e <= max_e; ++e) {
    BagScores tmpBags;
    score.emplace(e, tmpBags);
  }
  
  for (auto& ae : a->score) {
    for (auto& be : b->score) {
      uint32_t e = ae.first + be.first;
      if (e < min_e || e > max_e) continue;
      for (auto& a_bagmask : ae.second) {
        if (score.at(e).find(a_bagmask.first) == score.at(e).end()) {
          MaskScores tmpScores;
          score.at(e).emplace(a_bagmask.first, tmpScores);
        }
        if (be.second.find(a_bagmask.first) == be.second.end()) continue;
//        std::cerr << "Join jobs @" << min_e << "-" << e << "-" << max_e << ": " << be.second.at(a_bagmask.first).size() * a_bagmask.second.size() << std::endl;
        for (auto& a_maskScore : a_bagmask.second) {
          for (auto& b_maskScore : be.second.at(a_bagmask.first)) {
            Score s = a_maskScore.first + b_maskScore.first;
            if (s > parent->upperLimit_) continue;
            ForgetMask placeMask = a_maskScore.second | b_maskScore.second;
            if (opt_::MAXIMUM_RESULT_COUNT == 0
                || score.at(e).at(a_bagmask.first).count(s) < opt_::MAXIMUM_RESULT_COUNT) {
              score.at(e).at(a_bagmask.first).emplace(s, placeMask);
            }
          }
        }
      }
    }
  }
  if (opt_::FPT::MINIMUM_PROPAGATION_DEPTH > 0) {
    for (auto& e : score) {
      for (auto& bm : e.second) {
        std::set<Score> minScores;
        for (auto& ms : bm.second) {
          if (minScores.size() < opt_::FPT::MINIMUM_PROPAGATION_DEPTH) {
            minScores.emplace(ms.first);
          } else {
            Score maxScore = *std::max_element(minScores.begin(), minScores.end());
            if (ms.first < maxScore) minScores.emplace(ms.first);
            if (minScores.size() > opt_::FPT::MINIMUM_PROPAGATION_DEPTH)
              minScores.erase(maxScore);
          }
        }
        Score maxScore = *std::max_element(minScores.begin(), minScores.end());
        for (auto it = bm.second.begin(); it != bm.second.end();) {
          if (it->first > maxScore) it = bm.second.erase(it);
          else ++it;
        }
      }
    }
  }
  
//  if (opt_::MAXIMUM_RESULT_COUNT > 0) {
//    uint32_t x = opt_::MAXIMUM_RESULT_COUNT;
//    for (auto& e : score) {
//      for (auto& bm : e.second) {
//        std::multiset<Score> seenScores;
//        for (auto it = bm.second.begin(); it != bm.second.end();) {
//          seenScores.emplace(it->second);
//          if (seenScores.count(it->second) > x) {
//            it = bm.second.erase(it);
//          } else {
//            ++it;
//          }
//        }
//      }
//    }
//  }
  
  a->score.clear();
  b->score.clear();
}

ElnDist FPTOptimisation::VertMaskToElnDist(const VertMask& m) {
  ElnDist actualDist(parent_->possibleLocations_.size());
  for (size_t i = 0, j = 0; i < placedMask_.size(); ++i) {
    if (placedMask_[i] && m[i]) actualDist.set(j);
    if (placedMask_[i]) ++j;
  }
  return actualDist;
}

Score FPTOptimisation::ScoreVertex(MolVertPair v, VertMask f)
{  
  Score s = opt_::INF;
  const VertMask placed = f & placedMask_;
  if (placed.count() > parent_->electronsToAdd_) return s;
  
  ElnDist actualDist = VertMaskToElnDist(placed);
  
  ElnVertex u = mvp2ev_.at(v);
  SetElectronDistribution(actualDist);
  DetermineFormalCharges();
  s = CalculateVertexEnergy(u);
  return s;
}

String TDVertScore::VertMaskToNiceString(VertMask m) {
  std::ostringstream ss;
  ss << "[";
  VertMask placed = m; // & parent->actualPlacedMask_;
  MolecularGraph_p mg = parent->td_->GetSourceGraph()->GetSourceGraph()->GetSourceGraph();
  for (auto& mvp2vm : parent->pairMasks_) {
    if (!(placed & mvp2vm.second).count()) continue;
    uid_t a = mg->GetVertexIndex(mvp2vm.first.first);
    uid_t b = mg->GetVertexIndex(mvp2vm.first.second);
    if (a <= b) ss << "(" << a << "," << b;
    else ss << "(" << b << "," << a;
    ss << ";" << (placed & mvp2vm.second).count() << ") ";
  }
  ss << "]";
  return ss.str();
}

String TDVertScore::ToString() {
  std::ostringstream ss;
  for (auto& e : score) {
    ss << "Forgotten " << e.first << std::endl;
    for (auto& m : e.second) {
      ss << "BagState: " << VertMaskToNiceString(m.first - fMask);
//      ss << "  ForgotState: " << VertMaskToNiceString(m.first & fMask);
//      ss << "  Score: " << m.second << std::endl;
    }
  }
  return ss.str();
}

String TDVertScore::KindToString() {
  std::ostringstream ss;
  MolecularGraph_p mg = parent->td_->GetSourceGraph()->GetSourceGraph()->GetSourceGraph();
  if (tdProperties->kind.first == 'L') ss << "LEAF";
  if (tdProperties->kind.first == 'I') {
    ss << "INTRODUCE: ";
    uid_t a = mg->GetVertexIndex(tdProperties->kind.second.first);
    uid_t b = mg->GetVertexIndex(tdProperties->kind.second.second);
    if (a <= b) ss << a << "," << b;
    else ss << b << "," << a;
  }
  if (tdProperties->kind.first == 'F') {
    ss << "FORGET: ";
    uid_t a = mg->GetVertexIndex(tdProperties->kind.second.first);
    uid_t b = mg->GetVertexIndex(tdProperties->kind.second.second);
    if (a <= b) ss << a << "," << b;
    else ss << b << "," << a;
  }
  if (tdProperties->kind.first == 'J') {
    ss << "JOIN";
  }
  if (tdProperties->kind.first == 'R') {
    ss << "ROOT FORGET: ";
    uid_t a = mg->GetVertexIndex(tdProperties->kind.second.first);
    uid_t b = mg->GetVertexIndex(tdProperties->kind.second.second);
    if (a <= b) ss << a << "," << b;
    else ss << b << "," << a;
  }
  return ss.str();
}
