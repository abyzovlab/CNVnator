// C/C++ includes
#include <iostream>
#include <fstream>

#ifdef USE_YEPPP
#include <yepMath.h>
#include <yepLibrary.h>
#endif

using namespace std;

// Application includes
#include "AliParser.hh"
#include "HisMaker.hh"
#include "EXOnator.hh"
#include "IO.hh"
#include "Visualizer.hh"

int main(int argc,char *argv[])
{
  string usage = "\nCNVnator ";
#ifdef CNVNATOR_VERSION
  usage += CNVNATOR_VERSION;
#else
  usage += "v???";
#endif
  usage += "\n\nUsage:\n";
  usage += argv[0];
  usage += " -root out.root  [-genome name] [-chrom 1 2 ...] -tree  file1.bam ... [-lite]\n";
  usage += argv[0];
  usage += " -root out.root  [-genome name] [-chrom 1 2 ...] -merge file1.root ...\n";
  usage += argv[0];
  usage += " -root file.root  [-genome name] [-chrom 1 2 ...] -vcf [file.vcf.gz | file.vcf] [-rmchr] [-addchr]\n";
  usage += argv[0];
  usage += " -root file.root  [-genome name] [-chrom 1 2 ...] -idvar [file.vcf.gz | file.vcf] [-rmchr] [-addchr]\n";
  usage += argv[0];
  usage += " -root file.root  [-genome name] [-chrom 1 2 ...] -mask strict.mask.file.fa.gz [-rmchr] [-addchr]\n";
  usage += argv[0];
  usage += " -root file.root [-genome name] [-chrom 1 2 ...] [-d dir | -fasta file.fa.gz] -his bin_size\n";
  usage += argv[0];
  usage += " -root file.root [-genome name] [-chrom 1 2 ...] -baf bin_size [-hap] [-useid] [-nomask]\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -stat      bin_size\n";
  usage += argv[0];
  usage += " -root file.root                  -eval      bin_size\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -partition bin_size [-ngc]\n";
  // usage += argv[0];
  //usage += " -root file.root [-chrom 1 2 ...] -spartition bin_size [-gc]\n";
  usage += argv[0];
  usage += " -root file.root [-chrom 1 2 ...] -call      bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -root file.root -genotype bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -root file.root -view     bin_size [-ngc]\n";
  usage += argv[0];
  usage += " -pe   file1.bam ... -qual val(20) -over val(0.8) [-f file]\n";
  usage += argv[0];
  usage += "-root file.root [-chrom 1 2 ...] -cptrees newfile.root\n";
  usage += argv[0];
  usage += "-root file.root -ls\n";
  usage += "\n";
  usage += "Valid genomes (-genome option) are: NCBI36, hg18, GRCh37, hg19, mm9, hg38, GRCh38\n";

  if (argc < 2) {
    cerr<<"Not enough parameters."<<endl;
    cerr<<usage<<endl;
    return 0;
  }

#ifdef USE_YEPPP
  YepStatus yepStatus = yepLibrary_Init();
  if (yepStatus != YepStatusOk) {
    cerr<<"Yeppp library initialization failed with status "<<yepStatus<<"."<<endl;
    return 1;
  }
#endif

  static const int OPT_TREE       = 0x00001;
  static const int OPT_MERGE      = 0x00002;
  static const int OPT_HIS        = 0x00004;
  static const int OPT_HISMERGE   = 0x00008;
  static const int OPT_STAT       = 0x00010;
  static const int OPT_PARTITION  = 0x00020;
  static const int OPT_EPARTITION = 0x00040;
  static const int OPT_CALL       = 0x00080;
  static const int OPT_VIEW       = 0x00100;
  static const int OPT_VIEWER     = 0x00101;
  static const int OPT_GENOTYPE   = 0x00200;
  static const int OPT_EVAL       = 0x00400;
  static const int OPT_PE         = 0x00800;
  static const int OPT_PANEL      = 0x01000;
  static const int OPT_FIT        = 0x02000;

  static const int OPT_SPARTITION  = 0x04000;
  static const int OPT_PARTITION2D = 0x08000;
  static const int OPT_HIS_NEW     = 0x10000;
  static const int OPT_AGGREGATE   = 0x20000;
  
  static const int OPT_VCF = 0x50001;
  static const int OPT_IDVAR = 0x50002;
  static const int OPT_MASK = 0x50003;
  static const int OPT_BAF = 0x50004;
  static const int OPT_CALLBAF = 0x50005;

  static const int OPT_CPTREES = 0x60001;
  static const int OPT_LS = 0x60002;



  // tree, merge, his, stat, partition, spartition, call, view, genotype
  int max_opts = 10000, n_opts = 0, opts[max_opts], bins[max_opts], gbin = 0;
  for (int i = 0;i < n_opts;i++) bins[i] = 0;
  bool useGCcorr = true,useATcorr = false;
  bool lite = false,relaxCalling = false;
  string out_root_file(""),call_file(""),group_name(""),fastafile("");
  string chroms[1000],data_files[100000],root_files[100000] = {""},dir = ".";
  int n_chroms = 0,n_files = 0,n_root_files = 0,range = 128, qual = 20;
  double over = 0.8;
  double deltaAF = 0.25;
  Genome *genome = NULL;
  
  // vcf, idvar, mask option -rmchr -addchr
  bool rmchr=false,addchr=false;
  // baf option -hap -useid -nomask
  bool useHaplotype=false,useid=false,usemask=true;
  
  string signal="";
  string signal2="";


  int index = 1;
  while (index < argc) {
    string option = argv[index++];
    if (option == "-tree"  || option == "-merge" || option == "-pe" || option == "-vcf" || option == "-idvar" || option == "-mask" || option == "-cptrees") {
      if (option == "-tree")  opts[n_opts++] = OPT_TREE;
      if (option == "-merge") opts[n_opts++] = OPT_MERGE;
      if (option == "-pe")    opts[n_opts++] = OPT_PE;
      if (option == "-vcf")  opts[n_opts++] = OPT_VCF;
      if (option == "-idvar")    opts[n_opts++] = OPT_IDVAR;
      if (option == "-mask")    opts[n_opts++] = OPT_MASK;
      if (option == "-cptrees")    opts[n_opts++] = OPT_CPTREES;
      while (index < argc && argv[index][0] != '-')
	if (strlen(argv[index++]) > 0) data_files[n_files++] = argv[index - 1];
    } else if (option == "-his"        || option == "-his_new"     ||
	       option == "-hismerge"   ||
	       option == "-stat"       || option == "-eval"        ||
	       option == "-partition"  || option == "-partition2D" ||
	       option == "-epartition" || option == "-spartition"  ||
	       option == "-call"       || option == "-view"        || option == "-viewer"        ||
	       option == "-genotype"   || option == "-aggregate"   ||
         option == "-baf"        || option == "-callbaf") {
      int bs = 0;
      if (index < argc && argv[index][0] != '-') {
	TString tmp = argv[index++];
	if (!tmp.IsDigit()) {
	  cerr<<"Bin size must be integer for option '"<<option<<"'."<<endl;
	  cerr<<usage<<endl;
	  return 0;
	}
	bs = tmp.Atoi();
      }
      if (option == "-his")         opts[n_opts] = OPT_HIS;
      if (option == "-hismerge")    opts[n_opts] = OPT_HISMERGE;
      if (option == "-stat")        opts[n_opts] = OPT_STAT;
      if (option == "-partition")   opts[n_opts] = OPT_PARTITION;
      if (option == "-partition2D") opts[n_opts] = OPT_PARTITION2D;
      if (option == "-epartition")  opts[n_opts] = OPT_EPARTITION;
      if (option == "-spartition")  opts[n_opts] = OPT_SPARTITION;
      if (option == "-call")        opts[n_opts] = OPT_CALL;
      if (option == "-callbaf")        opts[n_opts] = OPT_CALLBAF;
      if (option == "-view")        opts[n_opts] = OPT_VIEW;
      if (option == "-viewer")        opts[n_opts] = OPT_VIEWER;
      if (option == "-genotype")    opts[n_opts] = OPT_GENOTYPE;
      if (option == "-his_new")     opts[n_opts] = OPT_HIS_NEW;
      if (option == "-eval")        opts[n_opts] = OPT_EVAL;
      if (option == "-aggregate")   opts[n_opts] = OPT_AGGREGATE;
      if (option == "-baf")         opts[n_opts] = OPT_BAF;
      bins[n_opts++] = bs;
    } else if (option == "-ls") {
      opts[n_opts++] = OPT_LS;
    } else if (option == "-panel") {
      opts[n_opts++] = OPT_PANEL;
    } else if (option == "-fit")   {
      opts[n_opts++] = OPT_FIT;
      if (index < argc && argv[index][0] != '-') group_name = argv[index++];
      else {
	cout<<"Please provide name of a sample group."<<endl;
	return 0;
      }
    } else if (option == "-outroot") {
      if (index < argc && argv[index][0] != '-') out_root_file = argv[index++];
      else {
	cout<<"Please provide new root-file name."<<endl;
	return 0;
      }
    } else if (option == "-root") {
      while (index < argc && argv[index][0] != '-')
	if (strlen(argv[index++]) > 0)
	  root_files[n_root_files++] = argv[index - 1];
      if (n_root_files == 0) {
	cerr<<"Please provide root-file name."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
    } else if (option == "-outroot") {
      if (index < argc && argv[index][0] != '-') out_root_file = argv[index++];
      else {
	cout<<"Please provide new root-file name."<<endl;
	return 0;
      }
    } else if (option == "-chrom") {
      while (index < argc && argv[index][0] != '-')
	chroms[n_chroms++] = argv[index++];
      if (n_chroms == 0) {
	cerr<<"Provide chromosome names."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
    } else if (option == "-ngc") {
      useGCcorr = false;
    } else if (option == "-at") {
      useATcorr = true;
    } else if (option == "-genome") {
      if (index < argc)	genome = Genome::get(argv[index++]);
    } else if (option == "-d") {
      if (index < argc && argv[index][0] != '-')
	dir = argv[index++];
      else cerr<<"No directory is given."<<endl;
    } else if (option == "-fasta") {
      if (index < argc && argv[index][0] != '-')
  fastafile = argv[index++];
      else cerr<<"Fasta file is not given."<<endl;
    }else if (option == "-qual") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No quality value is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      TString tmp = argv[index++];
      if (!tmp.IsDigit()) {
	cerr<<"Quality value must be integer."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      qual = tmp.Atoi();
    } else if (option == "-over") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No fraction of overlap is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      TString tmp = argv[index++];
      if (!tmp.IsFloat()) {
	cerr<<"Fraction of overlap must be number."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      over = tmp.Atof();
    } else if (option == "-f") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No file name is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      call_file = argv[index++];
    } else if (option == "-lite") {
      lite = true;
    } else if (option == "-rmchr") {
      rmchr=true;
    } else if (option == "-addchr") {
      addchr=true;
    } else if (option == "-useid") {
      useid=true;
    } else if (option == "-nomask") {
      usemask=false;
    } else if (option == "-hap") {
      useHaplotype=true;
    } else if (option == "-range") {
      range = atoi(argv[index++]);
    } else if (option == "-signal") {
      while((index < argc) && argv[index][0] != '-') {
        if(signal!="") signal+=" ";
        signal += argv[index++];
      }
    } else if (option == "-signal2") {
      while((index < argc) && argv[index][0] != '-') {
        if(signal2!="") signal2+=" ";
        signal2 += argv[index++];
      }
    } else if (option == "-relax") {
      relaxCalling = true;
    } else if (option == "-deltaAF") {
      if (index >= argc || argv[index][0] == '-') {
	cerr<<"No delta value is provided."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      TString tmp = argv[index++];
      if (!tmp.IsFloat()) {
	cerr<<"Quality value must be real number."<<endl;
	cerr<<usage<<endl;
	return 0;
      }
      deltaAF = tmp.Atof();
    } else if (option[0] == '-') {
      cerr<<"Unknown option '"<<option<<"'.\n"<<endl;
    }
  }

  if (out_root_file.length() <= 0) out_root_file = root_files[0];
  if (out_root_file.length() <= 0)
    cerr<<"WARNING: no name of root-file provided."<<endl;

  for (int o = 0;o < n_opts;o++) {
    int option = opts[o];
    int bin = bins[o]; if (bin <= 0) bin = gbin;
    if (option == OPT_TREE) { // tree
      HisMaker maker(out_root_file,genome);
      maker.setDataDir(dir);
      maker.produceTrees(chroms,n_chroms,data_files,n_files,lite);
    }
    if (option == OPT_VCF) { // vcf
      HisMaker maker(out_root_file,genome);
      maker.setDataDir(dir);
      maker.addVcf(chroms,n_chroms,data_files,n_files,rmchr,addchr);
    }
    if (option == OPT_IDVAR) { // idvar
      HisMaker maker(out_root_file,genome);
      maker.setDataDir(dir);
      maker.IdVar(chroms,n_chroms,data_files,n_files,rmchr,addchr);
    }
    if (option == OPT_MASK) { // mask
      HisMaker maker(out_root_file,genome);
      maker.setDataDir(dir);
      maker.MaskVar(chroms,n_chroms,data_files,n_files,rmchr,addchr);
    }
    if (option == OPT_CPTREES) { // cptrees
      IO io(out_root_file);
      if(n_files==1) io.cptrees(data_files[0],chroms,n_chroms);
      else cerr << "Provide one new root file name!" << endl;
    }
    if (option == OPT_MERGE) { // merge
      HisMaker maker(out_root_file,genome);
      maker.mergeTrees(chroms,n_chroms,data_files,n_files);
    }
    if (option == OPT_HIS ||
	option == OPT_HISMERGE) { // his
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      if(fastafile!="") {
        maker.setFastaFile(fastafile);
        maker.produceHistogramsFa(chroms,n_chroms,root_files,n_root_files,false);
        if (option == OPT_HISMERGE)  maker.produceHistogramsFa(chroms,n_chroms,root_files,n_root_files,true);
      } else {
        maker.setDataDir(dir);
        maker.produceHistograms(chroms,n_chroms,root_files,n_root_files,false);
        if (option == OPT_HISMERGE)  maker.produceHistograms(chroms,n_chroms,root_files,n_root_files,true);
      }
    }
    if (option == OPT_BAF) { // baf
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.produceBAF(chroms,n_chroms,useGCcorr,useHaplotype,useid,usemask);
    }
    if (option == OPT_STAT) { // stat
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.stat(chroms,n_chroms,useATcorr);
    }
    if (option == OPT_PARTITION) { // partition
      unsigned int flag=(usemask?FLAG_USEMASK:0)|(useid?FLAG_USEID:0)|(useHaplotype?FLAG_USEHAP:0)|(useGCcorr?FLAG_GC_CORR:0);
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      if(signal=="") maker.partition(chroms,n_chroms,false,useATcorr,useGCcorr,false,range);
      else maker.partitionSignal(bin,signal,flag,chroms,n_chroms,false,false,range);
    }
    if (option == OPT_PARTITION2D) { // partition
      unsigned int flag=(usemask?FLAG_USEMASK:0)|(useid?FLAG_USEID:0)|(useHaplotype?FLAG_USEHAP:0)|(useGCcorr?FLAG_GC_CORR:0);
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
//      maker.partition2D(chroms,n_chroms,false,useATcorr,useGCcorr,false,range);
      if(signal!="" && signal2!="") maker.partitionSignal2D(bin,signal,signal2,flag,chroms,n_chroms,range);
      else maker.partitionSignal2D(bin,"RD","SNP i1",flag,chroms,n_chroms,range);
    }
    if (option == OPT_EPARTITION) { // exome partition
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.partition(chroms,n_chroms,false,useATcorr,useGCcorr,true,range);
    }
    if (option == OPT_CALL) { // call
      unsigned int flag=(usemask?FLAG_USEMASK:0)|(useid?FLAG_USEID:0)|(useHaplotype?FLAG_USEHAP:0)|(useGCcorr?FLAG_GC_CORR:0);
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      if(signal=="") maker.callSVs(chroms,n_chroms,useATcorr,useGCcorr,deltaAF);
      else maker.callSVsSignal(bin,signal,flag,chroms,n_chroms,deltaAF);
    }
    if (option == OPT_CALLBAF) { // callbaf
      unsigned int flag=(usemask?FLAG_USEMASK:0)|(useid?FLAG_USEID:0)|(useHaplotype?FLAG_USEHAP:0)|(useGCcorr?FLAG_GC_CORR:0);
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.callBAF(chroms,n_chroms,useGCcorr,useHaplotype,useid,usemask);
    }
    if (option == OPT_VIEW) { // view
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      TApplication theApp("App",0,0);
      maker.view(root_files,n_root_files,useATcorr,useGCcorr);
      theApp.Run();
    }
    if (option == OPT_VIEWER) { // viewer
      unsigned int flags=(usemask?FLAG_USEMASK:0)|(useid?FLAG_USEID:0)|(useHaplotype?FLAG_USEHAP:0)|(useGCcorr?FLAG_GC_CORR:0);
      Visualizer vis(root_files,n_root_files,bin,flags);
      TApplication theApp("App",0,0);
      vis.prompt();
      theApp.Run();
    }
    if (option == OPT_GENOTYPE) { // genotype
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      TApplication theApp("App",0,0);
      maker.genotype(root_files,n_root_files,useATcorr,useGCcorr);
      theApp.Run();
    }
    if (option == OPT_EVAL) { // eval
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.eval(root_files,n_root_files,useATcorr,useGCcorr);
    }
    if (option == OPT_PE) { // pe
      HisMaker maker("null",genome);
      if (call_file.length() > 0)
	maker.pe_for_file(call_file,data_files,n_files,over,qual);
      else {
	TApplication theApp("App",0,0);
	maker.pe(data_files,n_files,over,qual);
	theApp.Run();
      }
    }
    if (option == OPT_SPARTITION) { // spartition
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.partition(chroms,n_chroms,true,useATcorr,useGCcorr,false,range);
    }
    if (option == OPT_HIS_NEW) { // his_new
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.produceHistogramsNew(chroms,n_chroms);
    }
    // EXOnator options
    if (option == OPT_PANEL) { // panel
      EXOnator exonator(out_root_file);
      exonator.makeTables();
    }
    if (option == OPT_FIT) { // fit
      EXOnator exonator(out_root_file);
      exonator.fit(group_name);
    }
    if (option == OPT_AGGREGATE) { // aggregate
      HisMaker maker(out_root_file,bin,useGCcorr,genome);
      maker.setDataDir(dir);
      maker.aggregate(root_files,n_root_files,chroms,n_chroms);
    }
    if (option == OPT_LS) { // aggregate
      IO io(out_root_file);
      io.ls();
    }
  }

  return 0;
}
