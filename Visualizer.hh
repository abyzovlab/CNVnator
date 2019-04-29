#ifndef __VISUALIZER_HH__
#define __VISUALIZER_HH__

#include "IO.hh"

class Visualizer
{
string *files;
IO **io;
int n_files;
int bin;
int flags;
TCanvas *canvas;

void generateViewBAF(string chrom,int start,int end);
void drawHistograms(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* raw,TH1 *his,TH1 *hisc,TH1 *hisp,TH1 *hism);
void drawHistogram2D(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* his2d);
void drawHistogramsBAF(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* hiss,TH1 *hisc);

public:
  Visualizer(string *_files,int _n_files,int _bin,int _flags);
  ~Visualizer();
  
  void prompt();
};

#endif
