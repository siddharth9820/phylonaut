#include "SimpleCladeExtractor.hpp"
#include <newick.hpp>
#include <vector>
#include <unordered_set>

unordered_set<Clade> SimpleCladeExtractor::extract(TaxonSet& ts) {
  string line;

  vector<Bipartition > bipartitions;
  unordered_set<Clade> clades;

  
  while(!gtfile.eof()) {
    getline(gtfile, line);
    newick_to_clades(line, ts, clades);
  }
  Clade all_taxa(ts, ts.taxa_bs);
  clades.insert(all_taxa);
  for (Taxon t : ts) {
    Clade c(ts, t);
    clades.insert(c);
  }
  vector<Clade> complements;
  for (const Clade& c : clades) {
    if (clades.count(c.complement()) == 0) {
      complements.push_back(c.complement());
    }
  }
  for (const Clade& c : complements) {
    clades.insert(c);
  }
  return clades;
  
}
