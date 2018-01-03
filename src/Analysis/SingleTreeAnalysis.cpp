#include "SingleTreeAnalysis.hpp"

#include "../TripartitionScorer/TripartitionScorer.hpp"
#include <sstream>

stringstream SingleTreeAnalysis::toTree(TaxonSet& ts, clade_bitset& bs, vector<Clade>& clades, TripartitionScorer& tps) {
  auto& subclades = tps.get_subclades(bs)[0];
  stringstream ss;
  if (bs.popcount() == 2) {
    vector<Taxon> tv;
    for (Taxon t : bs) {
      tv.push_back(t);
    }
    ss << "(" << ts[tv[0]] << "," << ts[tv[1]] << ")";
  }
  
  else if (bs.popcount() == 1) {
    ss << ts[*(bs.begin())];
  }
  
  else {
    ss << "(" << toTree(ts, subclades.first, clades, tps).str() << "," << toTree(ts, subclades.second, clades, tps).str() << ")";
  }
  
  return ss;
}

string SingleTreeAnalysis::analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) {
  stringstream ss;
  ss << toTree(ts, ts.taxa_bs, clades, tps).str();
  ss << ";";
  return ss.str();
}

