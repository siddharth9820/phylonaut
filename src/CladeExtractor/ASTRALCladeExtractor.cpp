#include "ASTRALCladeExtractor.hpp"
#include <util/Logger.hpp>
#include <newick.hpp>
#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>


string remapped_treefile(istream& input, TaxonSet& ts) {
  string line;
  //  ifstream ifile(input);
  
  char name[] = "/tmp/fileXXXXXX";
  mkstemp(name);
  
  ofstream of(name);
  while(!input.eof()) {
    getline(input, line);
    if (line.size())
      of << map_newick_names(line, ts) << endl;
  }
  
  return name;
}


unordered_set<Clade> ASTRALCladeExtractor::extract(TaxonSet& ts) {
  stringstream clade_stream;
  string s = "java -jar " + astralpath + " -i " + remapped_treefile(gtfile, ts) + " -k searchspace_norun -o /dev/null";  
  
  if (exact) 
    s += " -x ";  
  if (limited) 
    s += " -p0 ";
  if (extragtfile.good())
    s += " -e " + remapped_treefile(extragtfile, ts);  
  // if (!Logger::isEnabled("DEBUG") )
  //   s += " 2> /dev/null";


  DEBUG << s << endl;
  
  char buffer[128];
  
  FILE* stream = popen(s.c_str(), "r");
  
  stringstream result;
  while(!feof(stream)) {
    if(fgets(buffer,128,stream) != NULL) {
      result << buffer;
    }
  }

  
  stringstream cladestream_mapped(result.str());
  string line;
  
  stringstream unmapped;

  unordered_set<Clade> clades;
  
  while (!cladestream_mapped.eof()) {
    getline(cladestream_mapped, line);
    string s = unmap_clade_names(line, ts);
    clades.insert(Clade(ts, s));
  }

  
  
  return clades;
}
