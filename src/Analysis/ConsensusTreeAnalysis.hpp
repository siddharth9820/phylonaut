#ifndef __CONSENSUSTREEANALYSIS_HPP__
#define __CONSENSUSTREEANALYSIS_HPP__

#include "Analysis.hpp"

class ConsensusTreeAnalysis : public Analysis {
  double level;
public:
  ConsensusTreeAnalysis(double level_) : level(level_) {}
  virtual string analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps);
  virtual bool   requireAllBest() { return true; };
  vector<Clade> compatible_clades(vector<Clade>& clades);
  vector<Clade> children(Clade& c, vector<Clade>& clades);
};

#endif
