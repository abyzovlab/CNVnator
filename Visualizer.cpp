#include "Visualizer.hh"

// C/C++ includes
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;

// ROOT includes
#include <TFrame.h>
#include <TKey.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TTimer.h>
#include <Getline.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TMath.h>
#include <TF1.h>
#include <TLine.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <Math/DistFunc.h>
#include <TPRegexp.h>
#include <THashTable.h>
#include <TGraph.h>

#include <readline/readline.h>
#include <readline/history.h>

const vector<string> Visualizer::vocabulary=vector<string>({"set RD", "set RD raw", "set RD partition", "set RD call", "set SNP", "set SNP likelihood", "set SNP likelihood call", "set SNP count", "show", "ls"});

Visualizer::Visualizer(string *_files,int _n_files,int _bin,int _flags) : files(_files),n_files(_n_files),bin(_bin),flags(_flags),main_canvas(NULL) {
  io=new IO*[n_files];
  for(int i=0;i<n_files;i++) io[i]=new IO(files[i]);
  signals["RD"]=true;
  signals["RD raw"]=true;
  signals["RD partition"]=true;
  signals["RD call"]=true;
  
  signals["SNP"]=false;
  
  signals["SNP likelihood"]=false;
  signals["SNP likelihood call"]=false;
  
  signals["SNP count"]=false;
  
  chroms=vector<string>({"1","15","16","18"});
}

Visualizer::~Visualizer() {
  
}

void Visualizer::print_options() {
  cout << endl;
  cout << "Panel 1" << endl;
  cout << "-------" << endl;
  cout << "  RD " << signals["RD"] << endl;
  cout << "  RD raw " << signals["RD raw"] << endl;
  cout << "  RD partition " << signals["RD partition"] << endl;
  cout << "  RD call " << signals["RD call"] << endl << endl;

  cout << "Panel 2" << endl;
  cout << "-------" << endl;
  cout << "  SNP " << signals["SNP"] << endl << endl;

  cout << "Panel 3" << endl;
  cout << "-------" << endl;
  cout << "  SNP likelihood " << signals["SNP likelihood"] << endl;
  cout << "  SNP likelihood call " << signals["SNP likelihood call"] << endl << endl;

  cout << "Panel 4" << endl;
  cout << "-------" << endl;
  cout << "  SNP count " << signals["SNP count"] << endl << endl;
}

void Visualizer::setOption(string signal) {
  signals[signal]=true;
  cout << signal << " 1" << endl;
}

void Visualizer::unsetOption(string signal) {
  signals[signal]=false;
  cout << signal << " 0" << endl;
}

bool Visualizer::panel(int i) {
  if(i==1) return signals["RD"]||signals["RD row"]||signals["RD partition"]||signals["RD call"];
  if(i==2) return signals["SNP"];
  if(i==3) return signals["SNP likelihood"]||signals["SNP likelihood call"];
  if(i==4) return signals["SNP count"];
  return false;
}

int Visualizer::panels() {
  int ret=0;
  for(int i=1;i<5;i++) if(panel(i)) ret++;
  return ret;
}

int Visualizer::postoi(TString s) {
  if (s.IsDigit()) return s.Atoi();
  TString sp = s(0,s.Length() - 1);
  if (!sp.IsDigit()) return -1;
  char z = s[s.Length() - 1];
  if (z == 'M' || z == 'm') sp += "000000";
  if (z == 'k' || z == 'K') sp += "000";
  return sp.Atoi();
}

bool Visualizer::parseRegionOption(TString input,TString &chrom,
        int &start,int &end,TString &option)
{
  string tmp="";
  istringstream sin(input.Data());
  sin>>tmp;
  sin>>option;
  if(tmp=="") return false;
  int pos=tmp.find(":");
  if(pos==string::npos) return false;
  chrom=tmp.substr(0,pos);
  tmp.erase(0,pos+1);
  pos=tmp.find("-");
  if(pos==string::npos) return false;
  start=postoi(tmp.substr(0,pos));
  if(start==-1) return false;
  tmp.erase(0,pos+1);
  pos=tmp.find(" ");
  end=postoi(tmp);
  if(end==-1) return false;
  return true;
}

bool Visualizer::parseCommand(TString input) {
  string tmp="",tmp1="",sig="";
  istringstream sin(input.Data());
  sin>>tmp;
  while(sin>>tmp1) sig+=tmp1+" ";
  sig.pop_back();
  if(tmp=="show") {
    print_options();
  } else if(tmp=="set") {
    setOption(sig);
    return true;
  } else if(tmp=="unset") {
    unsetOption(sig);
    return true;
  } else if(tmp=="ls") {
    for(int i=0;i<n_files;i++) io[i]->ls();
  } else if(tmp=="gview") {
    gview();
  }
  return false;
  
}

char** Visualizer::completer(const char* text, int start, int end) {
  rl_attempted_completion_over = 1;
  return rl_completion_matches(text, Visualizer::completion_generator);
}

char* Visualizer::completion_generator(const char* text, int state) {
//  std::vector<std::string> vocabulary{"set RD", "set RD raw", "set RD partition", "set RD call", "set SNP", "set SNP likelihood", "set SNP likelihood call", "set SNP count", "show"};
  static std::vector<std::string> matches;
  static size_t match_index = 0;
  if (state == 0) {
    matches.clear();
    match_index = 0;
    std::string textstr = std::string(text);
    for (auto word : vocabulary) {
      if (word.size() >= textstr.size() &&
          word.compare(0, textstr.size(), textstr) == 0) {
        matches.push_back(word);
      }
    }
  }
  if (match_index >= matches.size()) {
    return nullptr;
  } else {
    return strdup(matches[match_index++].c_str());
  }
  return nullptr;
}

void Visualizer::prompt() {
  TTimer  *timer = new TTimer("gSystem->ProcessEvents();",50,kFALSE);
//  TString input = "";
  string input="";
  char* buf;
  rl_attempted_completion_function = completer;
  rl_basic_word_break_characters = (char*)"";
  while (input != "exit" && input != "quit") {
    TString chrom = "",option="";
    int start,end;
    if(input!="") if(!parseCommand(input)) if(parseRegionOption(input,chrom,start,end,option)) generateView(static_cast<string>(chrom),start,end);
    timer->TurnOn();
    timer->Reset();
    buf = readline("cnvnator> ");
    input=string(buf);
    if(buf!=nullptr) add_history(buf);
    free(buf);
//    input = Getline(">");
//    getline(cin,input);
//    input.ReplaceAll("\n","\0");
//    input = input.Remove(TString::kBoth,' ');
    timer->TurnOff();
  }
  delete timer;
}

void Visualizer::gview() {
//  int n=io[0]->nChromWithTree();
  int n=chroms.size();
  TCanvas *canvas = new TCanvas("canv","canv",1000,800);
  canvas->SetFillColor(kWhite);
  canvas->SetBorderMode(0); // No borders
  canvas->Divide(2,2);
  for(int i=0;i<n;i++) {
    TVirtualPad *pad = canvas->cd(i+1);
    cout<< "Chrom: " << chroms[i] << endl;
    generateView(chroms[i],1,io[0]->lenChromWithTree(chroms[i]),pad);
  }
  canvas->cd();
  canvas->Update();
}

void Visualizer::generateView(string chrom,int start,int end,TVirtualPad *canvas) {
  int win = 1*(end - start + 1);
  if (win < 10000) win = 10000;
  
//  if (canvas) delete canvas;
  TStyle *st = new TStyle("st","st");
  st->SetOptStat(false);  // No box with statistics
  st->SetOptTitle(false); // No box with title
  gROOT->SetStyle("st");
  if(!canvas) canvas = new TCanvas("canv","canv",1000,800);
  canvas->SetFillColor(kWhite);
  canvas->SetBorderMode(0); // No borders
  canvas->Divide(1,panels());
  
  TString title = chrom; title += ":";
  title += start; title += "-";
  title += end;
  canvas->SetTitle(title);
  int currp=0;
  
  // Panel 1
  if(panel(1)) {
    TH1 *raw=signals["RD raw"]?io[0]->getSignal(chrom,bin,"RD raw",0):NULL;
    TH1 *his=signals["RD raw"]?io[0]->getSignal(chrom,bin,"RD",0):NULL;
    TH1 *hisc=signals["RD"]?io[0]->getSignal(chrom,bin,"RD",flags):NULL;
    TH1 *hisp=signals["RD partition"]?io[0]->getSignal(chrom,bin,"RD partition",flags):NULL;
    TH1 *hism=signals["RD call"]?io[0]->getSignal(chrom,bin,"RD call",flags):NULL;
    TVirtualPad *pad = canvas->cd(++currp);
    pad->SetFillColor(kWhite);
    pad->SetLineColor(kWhite);
    pad->SetFrameLineColor(kWhite);
    pad->SetFrameBorderMode(0);
    drawHistograms(chrom,start,end,win,"RD",pad,raw,his,hisc,hisp,hism);
//  if (hisc) {
//    double mean,sigma;
//    TH1 *rd_his = getHistogram(getDistrName("chr1",bin_size,
//                                            useATcorr,useGCcorr));
//    getMeanSigma(rd_his,mean,sigma);
//    hisc->GetYaxis()->SetRangeUser(0,2*mean);
//    hisc->SetLineWidth(1);
//  }
  }

  if(panel(2)) {
    TH1 *his=io[0]->getSignal(chrom,bin,"SNP count",flags);
    TTree *tree=io[0]->getTree(chrom,"SNP");
    TVirtualPad *pad = canvas->cd(++currp);
    
    drawSNP(chrom,start,end,win,"SNP",pad,his,tree);
  }
  
  if(panel(3)) {
    TH1 *his=signals["SNP likelihood call"]?io[0]->getSignal(chrom,bin,"SNP likelihood call",flags):io[0]->getSignal(chrom,bin,"SNP likelihood",flags);
    TVirtualPad *pad = canvas->cd(++currp);
    pad->SetFillColor(kWhite);
    pad->SetLineColor(kWhite);
    pad->SetFrameLineColor(kWhite);
    pad->SetFrameBorderMode(0);
    drawHistogram2D(chrom,start,end,win,"likelihood",pad,his);
  }

  if(panel(4)) {
    TH1 *hisbafc=io[0]->getSignal(chrom,bin,"SNP count",flags);
    TVirtualPad *pad = canvas->cd(++currp);
    pad->SetFillColor(kWhite);
    pad->SetLineColor(kWhite);
    pad->SetFrameLineColor(kWhite);
    pad->SetFrameBorderMode(0);
    drawHistogramsBAF(chrom,start,end,win,"SNP count",pad,hisbafc,NULL);
  }

  canvas->cd();
  canvas->Update();
}

void Visualizer::drawHistograms(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* raw,TH1 *his,TH1 *hisc,TH1 *hisp,TH1 *hism)
{
  TH1 *main = hisc;
  if (!main) main = his;
  if (!main) main = raw;
  if (!main) return;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 0;
  for (int i = bs;i <= be;i++) {
    if (raw  && raw->GetBinContent(i)  > max) max = raw->GetBinContent(i);
    if (his  && his->GetBinContent(i)  > max) max = his->GetBinContent(i);
    if (hisc && hisc->GetBinContent(i) > max) max = hisc->GetBinContent(i);
  }
  max *= 1.05;

  main->Draw();
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);
  if (raw) {
    raw->Draw("same hist");
    raw->SetLineColor(kYellow);
  }
  if (his) {
    his->Draw("same hist");
    his->SetLineColor(kGray);
  }
  if (hisc) {
    hisc->Draw("same hist");
    hisc->SetLineColor(kBlack);
  }
  if (hisp) {
    hisp->Draw("same");
    hisp->SetLineColor(kRed);
    hisp->SetLineWidth(3);
  }
  if (hism) {
    hism->Draw("same hist");
    hism->SetLineColor(kGreen);
    hism->SetLineWidth(3);
  }
  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}

void Visualizer::drawHistogram2D(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* his2d)
{
  TH1 *main = his2d;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 1.05;

  main->Draw("CONT4");
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);

  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}

void Visualizer::drawHistogramsBAF(TString chrom,int start,int end,
            int win,TString title,
            TVirtualPad *pad,
            TH1* hiss,TH1 *hisc)
{
  TH1 *main = hiss;
  if (!main) main = hisc;

  int s  = start - win, e  = end + win;
  int n_bins = main->GetNbinsX();
  int bs = s/bin - 1; if (bs < 0) bs = 1;
  int be = e/bin + 1; if (be > n_bins) be = n_bins;
  double max = 0.5;

  main->Draw();
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  //main->GetYaxis()->SetRangeUser(0,max);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);
  if (hisc) {
    hisc->Draw("same hist");
    hisc->SetLineColor(kGreen);
    hisc->SetLineWidth(3);
  }
  TLine *line1 = new TLine(0,0,0,0),*line2 = new TLine(0,0,0,0);
  line1->SetX1(start); line1->SetX2(start);
  line1->SetY1(0);     line1->SetY2(max);
  line2->SetX1(end);   line2->SetX2(end);
  line2->SetY1(0);     line2->SetY2(max);
  line1->SetLineColor(kCyan);
  line2->SetLineColor(kCyan);
  line1->SetLineWidth(3);
  line2->SetLineWidth(3);
  line1->Draw(); line2->Draw();
}

void Visualizer::drawSNP(TString chrom,int start,int end,
                              int win,TString title,
                              TVirtualPad *pad,
                              TH1 *main,TTree *vcftree)
{
  int s  = start - win, e  = end + win;

  main->Draw();
  main->GetXaxis()->SetRangeUser(s,e);
  main->GetXaxis()->SetTitle(chrom);
  main->GetYaxis()->SetRangeUser(0,1.05);
  main->GetYaxis()->SetTitle(title);
  main->SetLineWidth(3);
  main->SetLineWidth(0);

  vcftree->SetMarkerSize(0.2);
  vcftree->SetMarkerStyle(21);
  vcftree->SetMarkerColor(kBlue);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==1 && (flag&2)","same");
  vcftree->SetMarkerColor(kCyan);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==1 && !(flag&2)","same");
  vcftree->SetMarkerColor(kRed);
  vcftree->Draw("nref/(nref+nalt):position","_gt%4==2 && (flag&2)","same");
  vcftree->SetMarkerColor(kOrange);
  vcftree->Draw("nref/(nref+nalt):position","_gt%4==2 && !(flag&2)","same");
  vcftree->SetMarkerColor(kGreen);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==0 && (flag&2)","same");
  vcftree->SetMarkerColor(kYellow);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==0 && !(flag&2)","same");
  vcftree->SetMarkerColor(kBlack);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==3 && (flag&2)","same");
  vcftree->SetMarkerColor(kGray);
  vcftree->Draw("nalt/(nref+nalt):position","_gt%4==3 && !(flag&2)","same");
}
