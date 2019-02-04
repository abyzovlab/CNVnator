#ifndef __FASTAPARSER__
#define __FASTAPARSER__

#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <zlib.h>

using namespace std;

// Parse genome data from .fasta.gz file
class FastaParser {
private:
  gzFile inFile;
  string buff;
  string c_name;
  string c_data;
  int n_cinb;
  bool eof;
  static const int sbuffs;
  char *sbuff;
  
public:
  inline string getName() {return c_name;} // returuns the name of current chromosome
  inline string getData() {return c_data;} // returns the the content of current chormosome
  bool nextChromosome(); // go to next chromosome in fasta file (returns false if eof)
  FastaParser(string fileName); // constructor - takes fa.gz filename as argument
  ~FastaParser(); // destructor
  
};

#endif
