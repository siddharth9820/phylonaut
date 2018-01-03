#ifndef __SCOREANALYSIS_HPP__
#define __SCOREANALYSIS_HPP__

#include "Analysis.hpp"
#include <Clade.hpp>

class ScoreAnalysis : public Analysis {
public:
  virtual string analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) {
    stringstream ss;
    Clade c(ts, ts.taxa_bs);
    ss << tps.adjust_final_score(tps.get_score(c));
    return ss.str();
  }
};

#endif
