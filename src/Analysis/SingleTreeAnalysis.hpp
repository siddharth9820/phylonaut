#ifndef __SINGLETREEANALYSIS_HPP__
#define __SINGLETREEANALYSIS_HPP__

#include "Analysis.hpp"
#include <Clade.hpp>

class SingleTreeAnalysis : public Analysis {
public:
  virtual string analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps);
  stringstream toTree(TaxonSet& ts, clade_bitset& bs, vector<Clade>& clades, TripartitionScorer& tps);
};

#endif
