#ifndef WASTRAL_HPP__
#define WASTRAL_HPP__

#include <vector>
#include <string>
#include "Clade.hpp"

class TripartitionScorer;
class TaxonSetExtractor;
class CladeExtractor;
class Analysis;

using namespace std;
class wASTRAL_ ;


struct Config {

  TripartitionScorer*  scorer;
  TaxonSetExtractor* taxon_extractor;

  vector<CladeExtractor*> extractors;  
  
  vector<Analysis*>    analyses;
  
  bool matrix;      //save the entire |X|^2 matrix
  
  string profile;     

  template<class InputIterator>
  void add_clades(InputIterator first, InputIterator last);
  vector<Clade>& get_clades();

  wASTRAL_* wASTRAL;
};

vector<string> wASTRAL(Config& conf);

class wASTRAL_ {
public:
  vector<string> run(Config& conf);  
  void uniqify_clades();  
  template<class InputIterator>
  void add_clades(InputIterator first, InputIterator last) {
    while (first !=last) {
      clades.push_back(*first);
      first++;
    }

    uniqify_clades();
    sort(clades.begin(), clades.end(), [](const Clade& a, const Clade& b){ return a.size() < b.size(); });

    

  }
  vector<Clade>& get_clades() { return clades; }

private:
  vector<Clade> clades;

    
};

template<class InputIterator>
void Config::add_clades(InputIterator first, InputIterator last)
{
  wASTRAL->add_clades(first, last);
}
#endif
