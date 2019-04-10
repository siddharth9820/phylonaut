#include "CountTreesAnalysis.hpp"
#include <vector>
#include <unordered_map>
#include <sstream>
#include "../TripartitionScorer/TripartitionScorer.hpp"
using namespace std;


unordered_map<clade_bitset, count_type> CountTreesAnalysis::run(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) {
  
  sort(clades.begin(), clades.end(), [](const Clade& a, const Clade& b){ return a.size() < b.size(); });

  unordered_map<clade_bitset, count_type> counts(clades.size());

  for (Clade& c : clades) {
    vector<pair<clade_bitset, clade_bitset> >& subclades_list = tps.get_subclades(c.get_taxa());

    if (c.size() <= 2 ) {
      counts[c.get_taxa()] = 1;
      continue;
    }
    
    count_type count = 0;
  
    for (auto& subclades : subclades_list) {      
      Clade c1(ts, subclades.first);
      Clade c2(ts, subclades.second);
      count += counts.at(c1.get_taxa()) * counts.at(c2.get_taxa());
    }
    
    counts[c.get_taxa()] = count;
   
    DEBUG << c.size() << "\t" << (double)counts[c.get_taxa()] << endl;
  }
  
  return counts;
}


string CountTreesAnalysis::analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) {
  unordered_map<clade_bitset, count_type> counts = run(ts, clades, tps);
  
  stringstream ss;
  DEBUG << ts.size() << "\t"  << (double)counts[ts.taxa_bs] << endl;
  ss << (double)counts[ts.taxa_bs]/(2 * ts.size() - 3);
  return ss.str();
}
