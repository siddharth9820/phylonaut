#ifndef TRIPARTITION_SCORER_HPP__
#define TRIPARTITION_SCORER_HPP__


#include <unordered_map>
#include <map>
#include <sstream>
#include <Clade.hpp>
#include <util/Logger.hpp>
#include <Quartet.hpp>

struct Config;

using namespace std;


typedef boost::multi_array<double, 2> twod_mat;

namespace std {
  template <typename T, typename U> struct hash<pair<T,U> > {
    size_t operator()(pair<T, U> x) {
      return hash<T>()(x.first) ^ hash<U>()(x.second);
    }
  };
};


class TripartitionScorer {
private:
  const TaxonSet* ts_;
  
  //Virtual methods
  //Do all the preprocessing
  virtual void setup(Config& conf, vector<Clade>& clades) {};
  
  //Returns the score of a tripartition
  virtual double score(const Tripartition& t)=0;


  //Return true if newscore is better than oldscore; by default
  // returns true if newscore < oldscore (i.e. minimizes the score)
  virtual bool better(double newscore, double oldscore);

public:
  
  //If the score should be adjusted before outputting; by default returns score
  virtual double adjust_final_score(double score);
  
  //Useful methods
  double get_score(Clade& clade);
  
  //used internally
  void init(Config& conf);
  
  vector<pair<clade_bitset, clade_bitset> >& get_subclades(const Clade& clade);
  vector<pair<clade_bitset, clade_bitset> >& get_subclades(const clade_bitset& clade);    

  TripartitionScorer () : ts_(NULL) {}
  
  void set_ts(const TaxonSet& newts) {
    DEBUG << clades.size() << endl;
    ts_ = &newts;
  }

  size_t clades_size() { return clades.size(); }
  
protected:
  const TaxonSet& ts() const {return *ts_;}

  
private:
  vector<Clade> clades;
  unordered_map <clade_bitset, size_t> clade_indices;
  double get_score(int clade_index);
  double get_score(Clade& clade, int clade_index);
  void set_score(size_t clade_index, double score, Clade& a1, Clade& a2); 
  
  vector <double> scores;
  vector <int> finished;
  vector <vector<pair<clade_bitset, clade_bitset> > > subclades;
  
  twod_mat* score_mat;
};






#endif
