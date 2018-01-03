#ifndef CLADE_EXTRACTOR__
#define CLADE_EXTRACTOR__
#include <vector>
#include <string>
#include <Clade.hpp>
#include <util/Options.hpp>

using namespace std;

class TaxonSetExtractor {
public:
  virtual TaxonSet& extract() = 0;  
};

class CladeExtractor {
public:
  virtual unordered_set<Clade> extract(TaxonSet& ts) = 0;  
};



#endif
