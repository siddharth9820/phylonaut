#ifndef SIMPLE_CLADE_EXTRACTOR__
#define SIMPLE_CLADE_EXTRACTOR__

#include "CladeExtractor.hpp"
#include <unordered_set>
#include <Clade.hpp>
#include <fstream>

using namespace std;

string findAstralJar();

class SimpleCladeExtractor : public CladeExtractor {
public:
  SimpleCladeExtractor(string& gtfile, bool root=true) :
    gtfile(*(new ifstream(gtfile))),
    root(root)
  {}

  SimpleCladeExtractor(istream& gtfile, bool root=true) :
    gtfile(gtfile),
    root(root)
  {}

  virtual unordered_set<Clade> extract(TaxonSet& ts);

private:
  istream& gtfile;
  bool root;

};


#endif
