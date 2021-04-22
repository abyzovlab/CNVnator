#include "IO.hh"
#include <algorithm>
#include <set>

const vector<string> IO::Signal = vector<string>({
  "RD",
  "RD unique",
  "RD raw",
  "RD l1",
  "RD l2",
  "RD l3",
  "RD partition",
  "RD call",
  "GC",
  "SNP",
  "SNP count",
  "SNP baf",
  "SNP maf",
  "SNP likelihood",
  "SNP i1",
  "SNP i2",
  "SNP i3",
  "SNP i4",
  "SNP likelihood partition",
  "SNP maf partition",
  "SNP i1 partition",
  "SNP likelihood call",
  "SNP maf call",
  "SNP i1 call",
});

// Distribution use same signal name + " dist"


IO::IO(string _rootfile) : rootfile(_rootfile),file(rootfile,"Update") {
  if (file.IsZombie()) {
    cerr<<"Can't open file '"<<rootfile<<"'."<<endl;
    return;
  }
  TIterator *it = file.GetListOfKeys()->MakeIterator();
  while (TKey *key = (TKey*)it->Next()) {
    TObject *obj = key->ReadObj();
    if (obj->IsA() != TTree::Class()) {
      string name = obj->GetName();
      string binstr="bin_";
      if (name.compare(0,binstr.size(),binstr)==0)
      bin_dir.push_back(name.substr(binstr.size()));
    } else {
      string chrom = obj->GetName();
      string title = obj->GetTitle();
      if (chrom == "") continue;
      int len;
      try {len=stoi(title.substr(title.find(";") + 1));}
      catch(...) {len=0;}
      
      string vcfstr="snp_";
      string vcfstr_old="vcf_";
      if (chrom.compare(0,vcfstr.size(),vcfstr)==0) {
        snp_tree.push_back(chrom.substr(vcfstr.size()));
        snp_tree_len.push_back(len);
      } else {
        if (chrom.compare(0,vcfstr_old.size(),vcfstr_old)!=0) {
          rd_tree.push_back(chrom);
          rd_tree_len.push_back(len);
        }
      }
    }
  }
}

IO::~IO() {
  file.Close();
}

TString IO::treeName(string chr,string signal) {
  if(signal=="RD") return chr;
  else if(signal=="SNP") return "snp_"+chr;
  cerr << "Unknown type of tree!" << endl;
  return "";
}

TString IO::signalName(string chr, int bin, string signal, int unsigned flags) {
  TString ret("");
  if(signal=="RD") {
    ret = "his_rd_p_" + chr + "_";
    ret += to_string(bin);
    if (flags&FLAG_AT_CORR) ret += "_AT";
    if (flags&FLAG_GC_CORR) ret += "_GC";
    return ret;
  } else if(signal=="RD unique") {
    ret = "his_rd_u_" + chr + "_";
    ret += to_string(bin);
    return ret;
  } else if(signal=="RD raw") {
    return signalName(chr,bin,"RD",flags)+"_raw";
  } else if(signal=="RD l1") {
    return signalName(chr,bin,"RD",flags)+"_l1";
  } else if(signal=="RD l2") {
    return signalName(chr,bin,"RD",flags)+"_l2";
  } else if(signal=="RD l3") {
    return signalName(chr,bin,"RD",flags)+"_l3";
  } else if(signal=="RD partition") {
    ret = signalName(chr,bin,"RD",0) + "_partition";
    if (flags&FLAG_AT_CORR) ret += "_AT";
    if (flags&FLAG_GC_CORR) ret += "_GC";
    return ret;
  } else if(signal=="RD call") {
    return signalName(chr,bin,"RD partition",flags)+"_merge";
  } else if(signal=="GC") {
    ret = chr + "_gc_";
    ret += bin;
    return ret;
  } else if(signal=="RD dist") {
    ret = "rd_p_";
    if (flags&FLAG_SEX) ret += "xy_";
    if (flags&FLAG_AT_CORR)  ret += "AT_";
    if (flags&FLAG_AT_CORR)  ret += "GC_";
    ret += bin;
    return ret;
  } else if(signal=="SNP") {
    return signalName(chr,bin,"SNP count",flags);
  } else if(signal=="SNP count") {
    ret = "snp_bafc_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP baf") {
    ret = "snp_baf_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP maf") {
    ret = "snp_maf_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP maf partition") {
    ret = "snp_maf_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_partition";
  } else if(signal=="SNP maf call") {
    ret = "snp_maf_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_call";
  } else if(signal=="SNP likelihood") {
    ret = "snp_likelihood_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP likelihood partition") {
    ret = "snp_likelihood_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_partition";
  } else if(signal=="SNP likelihood call") {
    ret = "snp_likelihood_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_call";
  } else if(signal=="SNP i1") {
    ret = "snp_i1_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP i1 partition") {
    ret = "snp_i1_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_partition";
  } else if(signal=="SNP i1 call") {
    ret = "snp_i1_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret+"_call";
  } else if(signal=="SNP i2") {
    ret = "snp_i2_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP i3") {
    ret = "snp_i3_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP i4") {
    ret = "snp_i4_"+chr+"_";
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    return ret;
  } else if(signal=="SNP count dist") {
    ret = "snp_dist_bafc_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP baf dist") {
    ret = "snp_dist_baf_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP maf dist") {
    ret = "snp_dist_maf_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP likelihood dist") {
    ret = "snp_dist_likelihood_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP i1 dist") {
    ret = "snp_dist_i1_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP i2 dist") {
    ret = "snp_dist_i2_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP i3 dist") {
    ret = "snp_dist_i3_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  } else if(signal=="SNP i4 dist") {
    ret = "snp_dist_i4_"+(chr==""?"":chr+"_");
    ret += to_string(bin);
    if(flags&FLAG_USEMASK) ret += "_mask";
    if(flags&FLAG_USEID) ret += "_id";
    if(flags&FLAG_USEHAP) ret += "_hap";
    if(flags&FLAG_SEX) ret += "_sex";
    return ret;
  }
  return "";
}

int IO::lenChromWithTree(string chr,bool vcf) {
  if(!vcf) { for(int i=0;i<rd_tree.size();i++)
    if(rd_tree[i]==chr) return rd_tree_len[i];
  } else { for(int i=0;i<snp_tree.size();i++)
    if(snp_tree[i]==chr) return snp_tree_len[i];
  }
  return -1;
}

int IO::nChromWithTree(bool vcf) {
  if(!vcf) return rd_tree.size();
  else return snp_tree.size();
}

int IO::getChromNamesWithTree(string *chrom_names,bool vcf) {
  if(!vcf) {
    for(int i=0;i<rd_tree.size();i++) chrom_names[i]=rd_tree[i];
    return rd_tree.size();
  }
  else {
    for(int i=0;i<snp_tree.size();i++) chrom_names[i]=snp_tree[i];
    return snp_tree.size();
  }
}

string IO::chromWithTree(int i,bool vcf) {
  if(!vcf) return rd_tree[i];
  else return snp_tree[i];
}

int IO::chromLenWithTree(int i,bool vcf) {
  if(!vcf) return rd_tree_len[i];
  else return snp_tree_len[i];
}

TString IO::getDirName(int bin) {
  TString ret = "bin_";
  ret += bin;
  return ret;
}

TTree* IO::getTree(string chr,string signal) {
  return (TTree*)file.Get(treeName(chr,signal));
}

TH1* IO::getSignal(string chr, int bin, string signal, unsigned int flags) {
  TDirectory *d = NULL;
  if (bin > 0) {
    TString dir=getDirName(bin);
    d = (TDirectory*)file.Get(dir);
    if (!d) {
      cerr<<"Can't find directory '"<<dir<<"'."<<endl;
      return NULL;
    }
    d->cd();
  } else d = &file;
  TH1 *his = (TH1*)d->Get(signalName(chr,bin,signal,flags));
  if (!his) return NULL;
  gROOT->cd();
  TH1 *ret = (TH1*)his->Clone(his->GetName());
  return ret;
}

TH1* IO::newSignalTH1(string title,string chr, int bin, string signal, unsigned int flags,int res,double min,double max) {
  TString h_name = signalName(chr,bin,signal,flags);
  TString h_title = title; if(chr!="") {h_title += " "; h_title +=chr;}
  TH1 *ret = new TH1D(h_name,h_title,res,min,max);
  ret->SetDirectory(0);
  return ret;
}

TH2* IO::newSignalTH2(string title,string chr, int bin, string signal, unsigned int flags,int res1,double min1,double max1,int res2,double min2,double max2) {
  TString h_name = signalName(chr,bin,signal,flags);
  TString h_title = title; if(chr!="") {h_title += " "; h_title +=chr;}
  TH2 *ret = new TH2D(h_name,h_title,res1,min1,max1,res2,min2,max2);
  ret->SetDirectory(0);
  return ret;
}

bool IO::writeHistogramsToBinDir(int bin,TH1 *his1,TH1 *his2,TH1 *his3,TH1 *his4,TH1 *his5,TH1 *his6) {
  //TFile file(rootfile,"Update");
  //if(file.IsZombie()) {
  //  cerr << "Can't open file '" << rootfile << "'." <<endl;
  //  return false;
  //}
  if(bin>0) {
    file.cd();
    TString dir_name=getDirName(bin);
    TDirectory *dir = (TDirectory*)file.Get(dir_name);
    if (!dir) {
      cout<<"Making directory "<<dir_name<<" ..."<<endl;
      dir = file.mkdir(dir_name);
      dir->Write(dir_name);
    }
    if (!dir) {
      cerr<<"Can't find/create directory '"<<dir_name<<"'."<<endl;
      return false;
    }
    dir->cd();
  }
  if (his1) his1->Write(his1->GetName(),TObject::kOverwrite);
  if (his2) his2->Write(his2->GetName(),TObject::kOverwrite);
  if (his3) his3->Write(his3->GetName(),TObject::kOverwrite);
  if (his4) his4->Write(his4->GetName(),TObject::kOverwrite);
  if (his5) his5->Write(his5->GetName(),TObject::kOverwrite);
  if (his6) his6->Write(his6->GetName(),TObject::kOverwrite);
  //file.Close();
  return true;
}

void IO::ls() {
  cout << endl << "Root file: " << rootfile << endl << "-------------------" << endl;
  if(rd_tree.size()>0) {
    cout << "RD TREES:";
    for(int i=0;i<rd_tree.size();i++)
    cout << ((i==0)?" ":", ") << rd_tree[i];
  }
  cout << endl << endl;
  if(snp_tree.size()>0) {
    cout << "SNP TREES:";
    for(int i=0;i<snp_tree.size();i++)
    cout << ((i==0)?" ":", ") << snp_tree[i];
  }
  cout << endl << endl;
  if(bin_dir.size()>0) {
    cout << "DIRECTORIES FOR BIN SIZES:";
    for(int i=0;i<bin_dir.size();i++)
    cout << ((i==0)?" ":", ") << bin_dir[i];
  }
  cout << endl << endl;
}

void IO::cptrees(string newrootfile,string *user_chroms,int n_chroms) {
  set<string> userchr;
  TFile *newfile=new TFile(newrootfile.c_str(),"Recreate");
  for(int i=0;i<n_chroms;i++) userchr.insert(user_chroms[i]);
  for(int i=0;i<rd_tree.size();i++) if((n_chroms==0)||(userchr.find(rd_tree[i])!=userchr.end())) {
    cout << "Copying RD tree for chromosome " << rd_tree[i] << "." << endl;
    TTree* t = (TTree*)file.Get(treeName(rd_tree[i],"RD"));
    newfile->cd();
    TTree* nt=t->CloneTree();
    nt->Write(0,TObject::kOverwrite);
  }
  for(int i=0;i<snp_tree.size();i++) if((n_chroms==0)||(userchr.find(snp_tree[i])!=userchr.end())) {
    cout << "Copying SNP tree for chromosome " << snp_tree[i] << "." << endl;
    TTree* t = (TTree*)file.Get(treeName(snp_tree[i],"SNP"));
    newfile->cd();
    TTree* nt=t->CloneTree();
    nt->Write(0,TObject::kOverwrite);
  }
  
  //Trees with obsolite "vcf_" prefix
  TIterator *it = file.GetListOfKeys()->MakeIterator();
  while (TKey *key = (TKey*)it->Next()) {
    TObject *obj = key->ReadObj();
    if (obj->IsA() == TTree::Class()) {
      string chrom = obj->GetName();
      string title = obj->GetTitle();
      if (chrom == "") continue;
      int len;
      try {len=stoi(title.substr(title.find(";") + 1));}
      catch(...) {len=0;}
      string vcfstr="vcf_";
      if (chrom.compare(0,vcfstr.size(),vcfstr)==0) if((n_chroms==0)||(userchr.find(chrom.substr(4))!=userchr.end())) {
        cout << "Copying SNP tree for chromosome " << chrom << " (obsolete name vcf_).";
        TTree* t = (TTree*)file.Get(chrom.c_str());
        newfile->cd();
        TTree* nt=t->CloneTree();
        nt->SetName(chrom.replace(0,4,"snp_").c_str());
        nt->Write(0,TObject::kOverwrite);
        cout << " Tree saved as " << chrom << "." << endl;
      }
    }
  }
  newfile->Close();
  delete newfile;
}
