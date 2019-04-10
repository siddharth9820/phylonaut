#include "GLOBExtractor.hpp"
#include <newick.hpp>
#include <vector>
#include <unordered_set>
#include <util/Logger.hpp>


Bipartition compatible(Bipartition& b1, Bipartition&  b2) {
  Clade b1_set = b1.a1 + b1.a2;
  Clade b2_set = b2.a1 + b2.a2;
  Bipartition b1_r = Bipartition(b1.a1.overlap(b2_set), b1.a2.overlap(b2_set));
  Bipartition b2_r = Bipartition(b2.a1.overlap(b1_set), b2.a2.overlap(b1_set));
  if (b1_r.a1.contains(b2_r.a1) || b2_r.a1.contains(b1_r.a1) || b1_r.a1.overlap_size(b2_r.a1) == 0) {
    if (!(b1_r.a1==b2_r.a1 || b1_r.a2==b2_r.a1 || b1_r.a1==b2_r.a2 || b1_r.a2==b2_r.a2) ) {
      Clade b1_right_fail = b2.a1 - b2_r.a1;
      Clade b1_left_fail = b2.a2 - b2_r.a2;
      return Bipartition(b1_left_fail, b1_right_fail);
    }
  }
  return Bipartition(Clade(b1.a1.ts()), Clade(b1.a1.ts()));

}


void complete_unordered(Bipartition& bp, vector<Bipartition>& bipartitions, unordered_set<Clade>& completed, TaxonSet& ts) {
    Clade missing_taxa = Clade(ts, ts.taxa_bs) - bp.a1 - bp.a2;
    vector<double> losses(ts.size());
    for (Bipartition& other : bipartitions) {
      Bipartition add = compatible(bp, other);

      for (Taxon tax : add.a1) {
        losses[tax] += 1;
      }

      for (Taxon tax : add.a2) {
        losses[tax] -= 1;
      }
    }

    Bipartition c_bp(bp.a1, bp.a2);


    for (Taxon tax : missing_taxa) {
      if (losses[tax] <= 0) {
        c_bp.a1.add(tax);
      } else {
        c_bp.a2.add(tax);
      }
    }
    DEBUG << c_bp.a1.size() + c_bp.a2.size() << endl;

    if (c_bp.a1.size())
      completed.insert(c_bp.a1);
    if (c_bp.a2.size())
      completed.insert(c_bp.a2);
}

unordered_set<Clade> GLOBExtractor::extract(TaxonSet& ts) {
  string line;

  vector<Bipartition > bipartitions;
  unordered_set<Clade> clades;

// a1 get the bipartitions
  while(!gtfile.eof()) {
    getline(gtfile, line);

    Clade tree_taxa = newick_to_taxa(line, ts);

    clades.clear();
    newick_to_clades(line, ts, clades);
    for (const Clade& c : clades) {
      bipartitions.push_back(Bipartition(c, tree_taxa.minus(c)));
    }
  }

  unordered_set<Clade> completed;

  for (Bipartition bp : bipartitions) {
    complete_unordered(bp, bipartitions, completed, ts);
  }
  INFO << completed.size() << " clades completed" << endl;

  for (Taxon tax : ts.taxa_bs) {
    Clade c(ts);
    c.add(tax);
    completed.insert(c);
    completed.insert(c.complement());
  }

  for (const Clade& bp : completed) {
    DEBUG << bp.str () << endl;
  }
  return completed;
}
