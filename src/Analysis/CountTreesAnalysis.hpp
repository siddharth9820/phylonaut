#ifndef __COUNTTREESANALYSIS_HPP__
#define __COUNTTREESANALYSIS_HPP__

#include "Analysis.hpp"

typedef __float128 count_type;

class CountTreesAnalysis : public Analysis {
public:
  unordered_map<clade_bitset, count_type> run(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps);
  virtual string analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps);
  virtual bool   requireAllBest() { return true; };
};

#endif
