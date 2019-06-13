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
  
  return clades;
  
}
