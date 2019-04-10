#include "ConsensusTreeAnalysis.hpp"
#include "CountTreesAnalysis.hpp"

#include <util/Logger.hpp>
#include <sstream>
#include <cassert>
#include <vector>
#include "../TripartitionScorer/TripartitionScorer.hpp"

vector<Clade> ConsensusTreeAnalysis::compatible_clades(vector<Clade>& clades) {
  vector<Clade> compat;
  for (Clade& c : clades) {
    bool good = true;
    
    for (Clade& other : compat) {
      if (!c.compatible(other)) {
	good = false;
	break;
      }      
    }
    
    if (good) 
      compat.push_back(c);    
  }
  return compat;
}


vector<Clade> ConsensusTreeAnalysis::children(Clade& c, vector<Clade>& clades) {
  Clade remaining(c);

  vector<Clade> clist;

  if (c.size() == 1) {
    clist.push_back(c);
    return clist;
  }
  
  for (Clade& other : clades) {
    if (remaining.contains(other) && !(c == other)) {
      clist.push_back(other);
      remaining = remaining.minus(other);
    }
    if (remaining.size() == 0) {
      DEBUG << "Clade " << c.str() << " has " << clist.size() << " children" << endl;
      return clist;
    }
  }

  
  ERR << "Clade missing some subclades" << endl;
  ERR << remaining.size() << "/" << c.size() << endl;
  assert(false);
}

void clades_in_optimal_tree(clade_bitset& base, TripartitionScorer& tps,  unordered_set<clade_bitset>& clades) {
  if (clades.count(base))
    return;
  clades.insert(base);
  vector<pair<clade_bitset, clade_bitset> >& subclades_list = tps.get_subclades(base);
  for (auto& subclades : subclades_list) {
    clades_in_optimal_tree(subclades.first, tps, clades);
    clades_in_optimal_tree(subclades.second, tps, clades);
  }
}


string ConsensusTreeAnalysis::analyze(TaxonSet& ts, vector<Clade>& clades, TripartitionScorer& tps) {
  CountTreesAnalysis counter;
  unordered_map<clade_bitset, count_type> counts = counter.run(ts, clades, tps);
  unordered_map<clade_bitset, count_type> appearances;
  unordered_set<clade_bitset> in_optimal;

  clades_in_optimal_tree(ts.taxa_bs, tps, in_optimal);

  DEBUG << in_optimal.size() << endl;
  
  count_type total_trees = counts[ts.taxa_bs]/(2 * ts.size() - 3);


  
  vector<Clade> good_clades;
  for (Clade& c : clades) {
    appearances[c.get_taxa()] = counts[c.get_taxa()] * counts[c.complement().get_taxa()];
    
    if (c.get_taxa() == ts.taxa_bs || (appearances[c.get_taxa()] > (level * total_trees) || (level==1.0 && appearances[c.get_taxa()] == (total_trees)) || (level<=0))) {
      DEBUG << (double)appearances[c.get_taxa()] << "\t" << (double)total_trees << endl;
      if (in_optimal.count(c.get_taxa()))
	  good_clades.push_back(c);
    }
  }  

 
  sort(good_clades.begin(), good_clades.end(), [ & ](const Clade& a, const Clade& b){ return appearances[a.get_taxa()] > appearances[b.get_taxa()]; });
  DEBUG << "Found " << good_clades.size() << " clades that appear in at least " << level << " trees" << endl;
  
  vector<Clade> compatible = compatible_clades(good_clades);

  sort(compatible.begin(), compatible.end(), [ & ](const Clade& a, const Clade& b){ return a.size() > b.size(); });  
  DEBUG << "Found " << compatible.size() << " compatible clades" << endl;
  
  unordered_map<clade_bitset, string> newicks;
  unordered_map<clade_bitset, vector<Clade> > kids;  

  for (Clade& c : compatible ) {
    kids[c.get_taxa()] = children(c, compatible);
    DEBUG << c.size() << "\t" << kids[c.get_taxa()].size() << endl;
  }

  sort(compatible.begin(), compatible.end(), [ & ](const Clade& a, const Clade& b){ return a.size() < b.size(); });    


  for (Clade& c : compatible ) {

    stringstream ss;

    if (c.size() == 1) {
      newicks[c.get_taxa()] = ts[*c.begin()];
      continue;
    }

    double support = (double)(appearances[c.get_taxa()]/total_trees);

    DEBUG << c.str() << "\t" << support << "\t" << kids[c.get_taxa()].size() << endl;
    
    if (c.size() == 2) {
      
      vector<Taxon> tv;
      for (Taxon t : c) {
	tv.push_back(t);
      }
      ss << "(" << ts[tv[0]] << "," << ts[tv[1]] << "):" << support;
      newicks[c.get_taxa()] = ss.str();
      continue;
    }
    
    ss << "(";
    bool started=false;
    for (auto& k : kids[c.get_taxa()]) {
      if (started)
	ss << ",";
      started=true;
      ss << newicks[k.get_taxa()];      
    }
    ss << ")";
    newicks[c.get_taxa()] = ss.str();
  }

  return newicks[ts.taxa_bs] + ";";
}

