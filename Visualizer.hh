#ifndef __VISUALIZER_HH__
#define __VISUALIZER_HH__

#include "IO.hh"
#include <map>
#include <vector>
#include <string>


class Visualizer
{
string *files;
IO **io;
int n_files;
int bin;
int flags;
TCanvas *main_canvas;
map<string,bool> signals;
vector<string> chroms;

static const vector<string> vocabulary;

  
void gview();

void generateView(string chrom,int start,int end,TVirtualPad *canvas=NULL);

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
void drawSNP(TString chrom,int start,int end,
                              int win,TString title,
                              TVirtualPad *pad,
                              TH1 *main,TTree *vcftree);
  
bool parseRegionOption(TString input,TString &chrom,
        int &start,int &end,TString &option);
bool parseCommand(TString input);

int postoi(TString s);
void print_options();
void setOption(string signal);
void unsetOption(string signal);
bool panel(int i);
int panels();

//readline related functions
static char** completer(const char* text, int start, int end);
static char* completion_generator(const char* text, int state);

public:
  Visualizer(string *_files,int _n_files,int _bin,int _flags);
  ~Visualizer();
  
  void prompt();
};

#endif
