#ifndef ANALYSIS_HPP__
#define ANALYSIS_HPP__

#include <vector>
#include <TaxonSet.hpp>
#include <Clade.hpp>

class TripartitionScorer;

class Analysis {
protected: 
public:  
  
  virtual string analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) = 0;

  
  virtual bool   requireMatrix() { return false; };
  virtual bool   requireAllBest() { return false; };
};

#endif
