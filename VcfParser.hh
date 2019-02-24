#ifndef __VCFPARSER__
#define __VCFPARSER__

#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include "vcf.h"
#include "vcfutils.h"

using namespace std;

class VcfParser {
private:

  htsFile *f;
  bcf_hdr_t *h;
  bcf1_t *rec;
  
  bool stdin;
  int nsamples;
  string      *cnames_;
  int         *clens_;
  int          n_chr_;
  int cmax;
  int nrec;
  
  istream *in;
  
  string c_chr;
  int c_pos;
  string c_id;
  string c_ref;
  string c_alt;
  double c_qual;
  string c_filter;
  string c_info;
  string c_format;
  string c_rec;
  string c_gt;
  int c_nref;
  int c_nalt;
  void parseHeader();
  
public:
  inline string getChromosome() {return c_chr;}
  int getChromosomeIndex();
  inline int getPosition() {return c_pos;}
  inline string getId() {return c_id;}
  inline string getRef() {return c_ref;}
  inline string getAlt() {return c_alt;}
  inline int getNRef() {return c_nref;}
  inline int getNAlt() {return c_nalt;}
  inline double getQual() {return c_qual;}
  inline string getFilter() {return c_filter;}
  inline string getInfo() {return c_info;}
  inline string getFormat() {return c_format;}
  inline string getRecord() {return c_rec;}
  inline string getGT() {return c_gt;}
  int    numChrom() { return n_chr_; }
  string chromName(int i) { return (i >= 0 && i < n_chr_) ? cnames_[i] : ""; }
  int    chromLen(int i)  { return (i >= 0 && i < n_chr_) ? clens_[i] : 0; }
  
  VcfParser(string fileName);
  ~VcfParser();
  
  bool parseRecord(bool idvar=false);
};

#endif
