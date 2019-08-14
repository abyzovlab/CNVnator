#ifndef __IO_HH__
#define __IO_HH__

#define FLAG_SEX 0x0001
#define FLAG_GC_CORR 0x0010
#define FLAG_AT_CORR 0x0020
#define FLAG_USEMASK 0x0100
#define FLAG_USEID 0x0200
#define FLAG_USEHAP 0x0400


// C/C++ includes
#include <iostream>
#include <string>
#include <vector>
#include <map>
using namespace std;

// ROOT includes
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>

class IO {
  const static vector<string> Signal;
  const static map<string,unsigned int> Flag;

  TString rootfile;
  vector<string> rd_tree;
  vector<string> snp_tree;
  vector<string> bin_dir;
  vector<int> rd_tree_len;
  vector<int> snp_tree_len;
  
public:
  IO(string _rootfile);
  ~IO();

  TFile file;
  // NAMES
  TString treeName(string chr,string signal);
  TString signalName(string chr, int bin, string signal, unsigned int flags);
  inline TString distributionName(int bin, string signal, unsigned int flags) {
    return signalName("",bin,signal,flags);
  }
  TString getDirName(int bin);

  // CHROMOSOMES WITH TREES
  // Use snp insted vcf
  int lenChromWithTree(string chr,bool vcf=false);
  int nChromWithTree(bool vcf=false);
  int getChromNamesWithTree(string *chrom_names,bool vcf=false);
  string chromWithTree(int i,bool vcf=false);
  int chromLenWithTree(int i,bool vcf=false);
  
  // READ FROM ROOT FILE
  TTree* getTree(string chr,string signal);
  TH1* getSignal(string chr, int bin, string signal, unsigned int flags);
  inline TH1* getDistribution(int bin, string signal, unsigned int flags) {
    return getSignal("",bin,signal,flags);
  }

  // CREATE HISTOGRAMS
  // Remove TH1 and TH2 and use His
  TH1* newSignalTH1(string title,string chr, int bin, string signal, unsigned int flags,int res,double min,double max);
  inline TH1* newDistributionTH1(string title, int bin, string signal, unsigned int flags,int res,double min,double max) {
    return newSignalTH1(title,"",bin,signal,flags,res,min,max);
  }
  TH2* newSignalTH2(string title,string chr, int bin, string signal, unsigned int flags,int res1,double min1,double max1,int res2,double min2,double max2);
  inline TH2* newDistributionTH2(string title, int bin, string signal, unsigned int flags,int res1,double min1,double max1,int res2,double min2,double max2) {
    return newSignalTH2(title,"",bin,signal,flags,res1,min1,max1,res2,min2,max2);
  }

  // WRITE TO ROOT FILE
  bool writeHistogramsToBinDir(int bin,TH1 *his1 = NULL,TH1 *his2 = NULL, TH1 *his3 = NULL,TH1 *his4 = NULL, TH1 *his5 = NULL,TH1 *his6 = NULL);
  inline bool writeHistograms(TH1 *his1 = NULL,TH1 *his2 = NULL,TH1 *his3 = NULL,TH1 *his4 = NULL, TH1 *his5 = NULL,TH1 *his6 = NULL) {
    return writeHistogramsToBinDir(0,his1,his2,his3,his4,his5,his6);
  }
  
  //MISC
  void ls();
  void cptrees(string newrootfile,string *user_chroms=NULL,int n_chroms=0);
};

#endif
