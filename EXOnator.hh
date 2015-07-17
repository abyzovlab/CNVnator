#ifndef __EXONATOR_HH__
#define __EXONATOR_HH__

// C/C++ includes
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <algorithm>
#include "TFile.h"
#include "TTree.h"

// ROOT includes
#include <TMath.h>
#include <TF1.h>
#include <TH1D.h>

using namespace std;

// genomic interval
class Interval {
public:
	string ref, name, intv_str;
	int start, end;

	Interval(string ref, long start, long end, string name);
	string toString() const;
	bool operator<(const Interval& that) const;
};

// sample
class Sample {
public:
	string sample, group, cov_file;

	Sample(string sample, string group, string cov_file);
	string toString() const;
};

// count data
class ReadCount {
public:
	long rcnt, ndup_rcnt;

	ReadCount(long rcnt, long ndup_rcnt);
};

class EXOnator
{
private:
  static const string FIT_DESCRIPTION;

private:
  string _root_file;

public:
  EXOnator(string root_file) : _root_file(root_file) {}

public:
  void makeTables(string bed_file = "", string conf_file = "");
  void makeTables2(string bed_file = "", string conf_file = "");
  void fit(string group_name,bool bimodal = false);
};


#endif /* __EXONATOR_HH__ */
