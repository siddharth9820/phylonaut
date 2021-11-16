#include "TripartitionScorer.hpp"
#include "../CladeExtractor/CladeExtractor.hpp"
#include "../wASTRAL.hpp"
#include <util/Logger.hpp>


#include <limits>
#include <fstream>
#include <cmath>
#include <omp.h>
#include <iostream>

bool TripartitionScorer::better(double newscore, double oldscore) {
  return newscore < oldscore;
}

double TripartitionScorer::adjust_final_score(double score) {
  return score; 
}

double TripartitionScorer::get_score(int clade_index) {
  return get_score(clades[clade_index], clade_index);
}

double TripartitionScorer::get_score(Clade& clade) {
  return get_score(clade, clade_indices[clade.get_taxa()]);
}

double TripartitionScorer::get_score(Clade& clade, int clade_index) {
  double value = (double)nan("");

  
  
  if(finished[clade_index]){
    return scores[clade_index];
  }
  DEBUG << "SCORING TRIPARTITION " << clade_index << endl;

  if (clade.size() == 1) {   
    Clade eclade (ts());
    set_score(clade_index, 0, clade, eclade);
    return 0;
  }     
  
  for (size_t i = 0; i < clades.size(); i++) {
    Clade& subclade = clades[i];
    
    if (subclade.size() >= clade.size() || subclade.size() == 0 || !clade.contains(subclade)  )
      continue;
    
    Tripartition tp(ts(), clade, subclade);
    
    if (clade_indices.count(tp.a1.get_taxa()) == 0 || clade_indices.count(tp.a2.get_taxa()) == 0)
      continue;
    
    double current = score(tp)						\
      + get_score(tp.a1, clade_indices[tp.a1.get_taxa()]) + get_score(tp.a2, clade_indices[tp.a2.get_taxa()]);

    if (current == value) {
      if (tp.a1.get_taxa().ffs() < tp.a2.get_taxa().ffs())
	subclades[clade_index].push_back(make_pair(tp.a1.get_taxa(), tp.a2.get_taxa()));
    }
    if (std::isnan(value) || better(current, value) ) {
      value = current;
      set_score(clade_index, value, tp.a1, tp.a2);
    }
    
    
  }
  DEBUG << clade.str() << "\t" << get_subclades(clade).size() << endl;
  return value;
}


void TripartitionScorer::set_score(size_t clade_index, double score, Clade& a1, Clade& a2) {
  subclades[clade_index].clear();
  scores[clade_index] = score;
  finished[clade_index] = 1;
  if (a1.get_taxa().ffs() < a2.get_taxa().ffs())
    subclades[clade_index].push_back(make_pair(a1.get_taxa(), a2.get_taxa()));
  
  if (score_mat) {
    (*score_mat)[clade_index][clade_indices[a1.get_taxa()]] = score;
    (*score_mat)[clade_index][clade_indices[a2.get_taxa()]] = score;
  }  

}

vector<pair<clade_bitset, clade_bitset>>& TripartitionScorer::get_subclades(const Clade& clade) {
  return get_subclades(clade.get_taxa());
}
  
vector<pair<clade_bitset, clade_bitset>>& TripartitionScorer::get_subclades(const clade_bitset& clade) {
  if(clade_indices.count(clade) == 0){
    Clade c(ts(), clade);
    ERR << c.str() << " doesn't have subclades!" << endl;
    assert(false);
  }

  size_t ix = clade_indices[clade];
  
  
  return subclades.at(ix);
  
}


void TripartitionScorer::init(Config& conf) {
  DEBUG << clades.size() << endl;
  setup(conf, conf.get_clades());

  DEBUG << clades.size() << endl;
  
  for (Clade& c : conf.get_clades())
    clades.push_back(c);
       
  
  for (size_t i = 0; i < clades.size(); i++) {
    clade_indices[clades[i].get_taxa()] = i;
  }

  if (conf.matrix) {
    DEBUG << "Making score matrix" << endl;
    score_mat = new twod_mat(boost::extents[clades.size()][clades.size()]);
   
    double start = omp_get_wtime(); 
    #pragma omp parallel for collapse(2)
    for (size_t i = 0; i < clades.size(); i++) {
      for (size_t j = 0; j < clades.size(); j++) {
	(*score_mat)[i][j] = (double)nan("");
      }
    }
    double end = omp_get_wtime();
    std::cout << "Time to make score_mat is " << end-start << " seconds with " << omp_get_max_threads() << " threads" << std::endl;
  } else{
    std::cout << "not building matrix" << std::endl;
    score_mat = 0;
  }
	
  scores.resize(clades.size());
  finished.resize(clades.size());

  for (size_t i = 0; i < clades.size(); i++) {
    finished[i] = 0;
  }
  
  subclades.resize(clades.size());

}


