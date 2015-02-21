#ifndef __GENOMES__
#define __GENOMES__

// C/C++ includes
#include <string>
#include <ctype.h>
#include <iostream>
using namespace std;

class Genome
{
private:
  Genome(string name);

private:
  static const int NGS = 2;
  static Genome genomes[NGS];
  string gname_,other_gname_;
  static const int MAX_N_CHROMS_ = 100;
  static const string CHRX, CHRY;
  string cnames_[MAX_N_CHROMS_];
  int    clens_[MAX_N_CHROMS_];
  int    n_chr_;

public:
  static const string CHRALL, CHRSEX;
  static Genome *get(string name);
  static string canonicalChromName(string name);
  static string extendedChromName(string name);
  static bool   isSexChrom(string name);

public:
  string name()     { return gname_; }
  int    numChrom() { return n_chr_; }
  string chromName(int i) { return (i >= 0 && i < n_chr_) ? cnames_[i] : ""; }
  int    chromLen(int i)  { return (i >= 0 && i < n_chr_) ? clens_[i] : 0; }
  int    getChromosomeIndex(string chr);
};

#endif
