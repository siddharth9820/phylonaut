#ifndef DEFAULT_TAXON_SET_EXTRACTOR__
#define DEFAULT_TAXON_SET_EXTRACTOR__

#include <vector>
#include <unordered_set>
#include <string>
#include <fstream>
#include <Clade.hpp>
#include <newick.hpp>
#include <util/Options.hpp>
#include <util/Logger.hpp>
#include "CladeExtractor.hpp"

using namespace std;

class DefaultTaxonSetExtractor : public TaxonSetExtractor {
public:
  string treefile;
  TaxonSet* ts;
  DefaultTaxonSetExtractor(const string& treefile) : treefile (treefile) {
  }
  virtual TaxonSet& extract() {
    unordered_set<string> taxa;
    ifstream in(treefile);

    string s;
    int i = 0;
    while(!in.eof()) {
      DEBUG << "Tree " << i << endl;
      getline(in, s);
      if (s.size() == 0)
	continue;
      newick_to_ts(s, taxa);
      DEBUG << taxa.size() << endl;
      i++;
    }

    ts = new TaxonSet(taxa.size());

    for (const string& s : taxa) {
      ts->add(s);
      DEBUG << s << endl;
    }
    
    return *ts;
  }
  
};
#endif
