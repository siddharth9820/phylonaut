#ifndef ASTRAL_CLADE_EXTRACTOR__
#define ASTRAL_CLADE_EXTRACTOR__

#include "CladeExtractor.hpp"
#include <unordered_set>
#include <Clade.hpp>
#include <fstream>

using namespace std;

string findAstralJar();

class ASTRALCladeExtractor : public CladeExtractor {
public:
  ASTRALCladeExtractor(string astralpath, string gtfile, string extragtfile, bool exact=false, bool limited=false) :
    astralpath(astralpath),
    gtfile(*(new ifstream(gtfile))),
    extragtfile(*(new ifstream(extragtfile))),
    exact(exact),
    limited(limited)
  {}

  ASTRALCladeExtractor(string astralpath, istream& gtfile, bool exact=false, bool limited=false) :
    astralpath(astralpath),
    gtfile(gtfile),
    extragtfile(*(new ifstream(""))),
    exact(exact),
    limited(limited){}

  ASTRALCladeExtractor(string gtfile, string extragtfile="", bool exact=false, bool limited=false) :
    astralpath(findAstralJar()),
    gtfile(*(new ifstream(gtfile))),
    extragtfile(*(new ifstream(extragtfile))),
    exact(exact),
    limited(limited)  {}

  ASTRALCladeExtractor(istream& gtfile, bool exact=false, bool limited=false) :
    astralpath(findAstralJar()),
    gtfile(gtfile),
    extragtfile(*(new ifstream(""))),
    exact(exact),
    limited(limited){}


  virtual unordered_set<Clade> extract(TaxonSet& ts);

private:
  string astralpath;
  istream& gtfile;
  istream& extragtfile;
  bool exact;
  bool limited;
};


#endif
