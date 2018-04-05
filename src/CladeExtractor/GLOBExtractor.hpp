#ifndef GLOB_CLADE_EXTRACTOR__
#define GLOB_CLADE_EXTRACTOR__

#include "CladeExtractor.hpp"
#include <unordered_set>
#include <Clade.hpp>
#include <fstream>
using namespace std;


class GLOBExtractor : public CladeExtractor{
private:
  istream& gtfile;

public:

  GLOBExtractor(ifstream& gtfile) :
  gtfile(gtfile)
  {  }

  GLOBExtractor(string gtfile) :
  gtfile(*(new ifstream(gtfile)))
  {  }


  virtual unordered_set<Clade> extract(TaxonSet& ts);

};

#endif
