// Application includes
#include "EXOnator.hh"

const string EXOnator::FIT_DESCRIPTION = "EXOnator_fit";

// Interval
Interval::Interval(string ref, long start, long end, string name) {
	this->ref = ref;
	this->name = name;
	this->start = start;
	this->end = end;
	ostringstream strs;
	strs << ref << ":" << start << "-" << end;
	this->intv_str = strs.str();
}

string Interval::toString() const {
	ostringstream strs;
	strs << this->ref << ":" << this->start << "-" << this->end << ":"
			<< this->name;
	return strs.str();
}

// compare intervals
bool Interval::operator<(const Interval& that) const {
	bool res = false;
	if (this->ref < that.ref) {
		res = true;
	} else if (this->ref == that.ref) {
		if (this->start < that.start) {
			res = true;
		}
	}
	return res;
}

// Sample
Sample::Sample(string sample, string group, string cov_file) {
	this->sample = sample;
	this->group = group;
	this->cov_file = cov_file;
}

string Sample::toString() const {
	ostringstream strs;
	strs << this->sample << "\t" << this->group << "\t" << this->cov_file;
	return strs.str();
}

// ReadCount
ReadCount::ReadCount(long rcnt, long ndup_rcnt) {
	this->rcnt = rcnt;
	this->ndup_rcnt = ndup_rcnt;
}

// build root file
void EXOnator::makeTables(string bed_file,string conf_file)
{
  string files[] = {"T_dup.txt","N_dup.txt"};
  char buffer[5000];
  string chrom,gene,sid;
  int start,end;
  double vals[1000];
  for (int i = 0;i < 2;i++) {
    TFile iroot(_root_file.c_str(),"update");
    TTree tree(files[i].substr(0,5).c_str(),files[i].substr(0,5).c_str());
    tree.Branch("chrom",&chrom);
    tree.Branch("start",&start,"start/I");
    tree.Branch("end",  &end,  "end/I");
    tree.Branch("gene", &gene);

    ifstream fin(files[i].c_str());
    fin.getline(buffer,5000);
    istringstream line0(buffer);
    line0>>sid>>sid>>sid>>sid;
    double *address = vals;
    while (line0>>sid)
      tree.Branch(sid.c_str(),address++,(sid + "/D").c_str());

    while (fin.good()) {
      fin.getline(buffer,5000);
      istringstream line(buffer);
      line>>chrom>>start>>end>>gene;
      address = vals;
      while (line>>(*address)) address++;
      tree.Fill();
    }
    fin.close();

    tree.Write("",TObject::kOverwrite);
    iroot.Close();
  }
}

void EXOnator::makeTables2(string bed_file,string conf_file)
{
  if (bed_file.length() == 0)
    bed_file = "/data5/beutler/m069635/gpanel/tumor/alexej/counts/S0674682_Covered_b37.bed";
  if (conf_file.length() == 0)
    conf_file = "/data5/beutler/m069635/gpanel/tumor/alexej/counts/config.txt";

  // read the bed file [ref start end gene_name]
  vector<Interval> bed;
  ifstream ibed(bed_file.c_str());
  string line, ref, gene;
  int start, end;
  while (getline(ibed, line)) {
    if (line[0] != '#' && line.length() != 0) {
      istringstream ss(line);
      ss >> ref >> start >> end >> gene;
      bed.push_back(Interval(ref, start, end, gene));
    }
  }
  ibed.close();
  
  // sort the intervals and make interval list
  sort(bed.begin(), bed.end());
  
  // read the config file [sample group file]
  map<string, vector<Sample> > conf;
  ifstream iconf(conf_file.c_str());
  string sample_str, group, cov_file;
  while (getline(iconf, line)) {
    if (line[0] != '#' && line.length() != 0) {
      istringstream ss(line);
      ss >> sample_str >> group >> cov_file;
      if (conf.find(group) == conf.end()) {
	conf[group] = vector<Sample>();
      }
      conf[group].push_back(Sample(sample_str, group, cov_file));
    }
  }
  iconf.close();
  
  // load data into root [ref start end rcnt ndup_rcnt]
  TFile iroot(_root_file.c_str(),"update");
  TTree tree_dup((group + "_dup").c_str(), (group + "_dup").c_str());
  TTree tree_ndup((group + "_ndup").c_str(), (group + "_ndup").c_str());

  string t_ref, t_name;
  int    t_start, t_end;
  tree_dup.Branch("chrom",&t_ref);
  tree_dup.Branch("start",&t_start,"start/I");
  tree_dup.Branch("end",  &t_end,  "end/I");
  tree_dup.Branch("gene", &t_name);
  
  tree_ndup.Branch("chrom",&t_ref);
  tree_ndup.Branch("start",&t_start,"start/I");
  tree_ndup.Branch("end",  &t_end,  "end/I");
  tree_ndup.Branch("gene", &t_name);
  
  for(map<string, vector<Sample> >::iterator mitr = conf.begin(); mitr != conf.end(); mitr++) {
    // make trees for this group
    string group = mitr->first;
    
    // add interval list to tree
    for (vector<Interval>::iterator vitr = bed.begin(); vitr != bed.end(); vitr++) {
      t_ref = vitr->ref.c_str();
      t_start = vitr->start;
      t_end = vitr->end;
      t_name = vitr->name.c_str();
      tree_dup.Fill();
      tree_ndup.Fill();
    }
    
    // add samples to tree
    for(vector<Sample>::iterator vitr = (mitr->second).begin(); vitr != (mitr->second).end(); vitr++) {
      // read the count file
      string ref;
      long start, end, rcnt, ndup_rcnt;
      long tot_rcnt = 0, tot_ndup_rcnt = 0;
      ifstream icnt((vitr->cov_file).c_str());
      map<string, ReadCount> cmap;
      while (getline(icnt, line)) {
	istringstream ss(line);
	ss >> ref >> start >> end >> rcnt >> ndup_rcnt;
	ostringstream key;
	key << ref << ":" << start << "-" << end;
	cmap.insert(pair<string, ReadCount>(key.str(), ReadCount(rcnt, ndup_rcnt)));
	tot_rcnt += rcnt;
	tot_ndup_rcnt += ndup_rcnt;
      }
      icnt.close();
      
      // add to root
      string sample = vitr->sample;
      double t_rcnt, t_ndup_rcnt;
      tree_dup.Branch(sample.c_str(), &t_rcnt,
		      (vitr->sample + "/D").c_str());
      tree_ndup.Branch(sample.c_str(), &t_ndup_rcnt,
		       (vitr->sample + "/D").c_str());
      for (vector<Interval>::iterator intv = bed.begin(); intv != bed.end(); intv++) {
	map<string, ReadCount>::iterator cnt = cmap.find(intv->intv_str);
	if (cnt != cmap.end()) {
	  t_rcnt = cnt->second.rcnt / (1.0 * tot_rcnt);
	  tree_dup.Fill();
	  t_ndup_rcnt = cnt->second.ndup_rcnt / (1.0 * tot_ndup_rcnt);
	  tree_ndup.Fill();
	}
      }
    }
  }

  tree_dup.Write();
  tree_ndup.Write();
  iroot.Close();
}

double gaus1(double *x,double *par)
{
  double inv = 1; if (par[2] != 0) inv = 1./par[2];
  double arg = (x[0] - par[1])*inv;
  return 0.4*inv*par[0]*TMath::Exp(-0.5*arg*arg);
}

double gaus2(double *x,double *par)
{
  double inv1 = 1; if (par[2] != 0) inv1 = 1./par[2];
  double arg1 = (x[0] - par[1])*inv1;
  double val1 = 0.4*inv1*par[0]*TMath::Exp(-0.5*arg1*arg1);

  double inv2 = 1; if (par[5] != 0) inv2 = 1./par[5];
  double arg2 = (x[0] - par[4])*inv2;
  double val2 = 0.4*inv2*par[3]*TMath::Exp(-0.5*arg2*arg2);

  return val1 + val2;
}

double getRMS(double *arr,int n,double mean,double rms = 0)
{
  if (rms == 0) rms = mean;
  double min = mean - 2.5*rms;
  double max = mean + 2.5*rms;
  double ret = 0;
  int count = 0;
  if (n <= 0) return ret;
  for (int i = 0;i < n;i++)
    if (arr[i] > min && arr[i] < max) {
      ret += (arr[i] - mean)*(arr[i] - mean);
      count++;
    }
  if (count > 0) return sqrt(ret/count);
  else           return 0;
}

double getKolmogorovProb(double *arr,int n,TF1 *func)
{
  double max = 0;
  double norm1 = 1./n;
  double norm2 = 1./func->Integral(func->GetXmin(),func->GetXmax());
  for (int i = 0;i < n;i++) {
    double diff = TMath::Abs((i + 1)*norm1 -
			     func->Integral(arr[i],func->GetXmax())*norm2);
    if (diff > max) max = diff;
  }

  return TMath::KolmogorovProb(max*TMath::Sqrt(n));
}

double gausIntegral(TF1 *func,double x)
{
  double mean    = func->GetParameter(1);
  double sigma   = func->GetParameter(2);
  if (sigma <= 0) sigma = 0.1;
  double delta = (x - mean)/sigma*0.707;
  double ret   = 0;
  if (delta > 0) ret = 0.5*(1 - TMath::Erf( delta));
  else           ret = 0.5*(1 + TMath::Erf(-delta));
  return ret;
}

double getKolmogorovProb1(double *arr,int n,TF1 *func)
{
  double max     = 0;
  double norm1   = 1./n;
  double mean    = func->GetParameter(1);
  double sigma   = func->GetParameter(2);
  if (sigma <= 0) sigma = 0.1;
  double sigma_o = 1/sigma;
  for (int i = 0;i < n;i++) {
    double x   = (arr[i] - mean)*sigma_o*0.707;
    double val = 0;
    if (x > 0) val = 0.5*(1 - TMath::Erf( x));
    else       val = 0.5*(1 + TMath::Erf(-x));
    double diff = TMath::Abs((i + 1)*norm1 - val);
    if (diff > max) max = diff;
  }

  return TMath::KolmogorovProb(max*TMath::Sqrt(n));
}

double getMedian(double *arr,int n)
{
  if (n%2 == 0) return 0.5*(arr[n/2 - 1] + arr[n/2]);
  return arr[n/2];
}

double getMedianSec(double *arr,int n)
{
  for (int i = 0;i < n - 1;i++) {
    double v3 = arr[i] - 0.01*(arr[i] - arr[i + 1]);
    double v1 = 0.5*v3, v2 = 0.75*v3;
    int count = 0;
    int b0 = 0, b1 = 0, b2 = 0, b3 = 0;
    for (int j = 0;j < n;j++) {
      double v = arr[j];
      if        (v < v1) {
	count++;
	b0++;
      } else if (v1 <= v && v < v2) {
	b1++;
      } else if (v2 <= v && v < v3) {
	count++;
	b2++;
      } else if (v3 <= v) {
	b3++;
      }
    }
    //cout<<v3<<" "<<count<<" "<<b0<<" "<<b1<<" "<<b2<<" "<<b3<<endl;
    if (b0 > 1 && b1 > 1 && b2 > 1 && b3 > 1 && count < n/2) return v3;
  }
  return 0;
}


void EXOnator::fit(string group_name,bool bimodal)
{
  TFile file(_root_file.c_str(),"Update");
  if (file.IsZombie()) { 
    cerr<<"Can't open file '"<<_root_file<<"'."<<endl;
    return;
  }
  TTree *tree = (TTree*)file.Get(group_name.c_str());
  if (!tree) {
    cerr<<"Can't find tree for group '"<<group_name<<"' in file '"
	<<_root_file<<"'."<<endl;
    return;
  }

  int Ndesc   = 4, Nbranch = 0;
  TObjArray *brs = tree->GetListOfBranches();
  if (brs) Nbranch = brs->GetEntries();
  if (Nbranch < Ndesc) {
    cerr<<"Tree '"<<group_name<<"' is of wrong format."<<endl;
    return;
  }

  string *chrom = new string(""), *tmp = new string("");
  int N = Nbranch - Ndesc, *ind = new int[N], start, end;
  tree->SetBranchAddress(brs->At(0)->GetName(),&chrom);
  tree->SetBranchAddress(brs->At(1)->GetName(),&start);
  tree->SetBranchAddress(brs->At(2)->GetName(),&end);
  tree->SetBranchAddress(brs->At(3)->GetName(),&tmp);

  double *vals = new double[N], *arr = new double[N], *addr = vals;
  for (int i = Ndesc;i < Nbranch;i++) {
    if (brs->At(i)->GetName() == FIT_DESCRIPTION) break;
    tree->SetBranchAddress(brs->At(i)->GetName(),addr++);
  }

  double pars[30] =
    {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  TBranch *br = tree->GetBranch(FIT_DESCRIPTION.c_str());
  if (br != NULL) tree->GetListOfBranches()->Remove(br);
  br = tree->Branch(FIT_DESCRIPTION.c_str(),&pars,"pars[30]/D");

  int nb = 100;
  TH1 *his = new TH1D("his","his",nb,0,1);
  int n_ent = tree->GetEntries();
  cout<<"n_ent = "<<n_ent<<endl;
  for (int e = 0;e < n_ent;e++) {
    tree->GetEntry(e);

    ostringstream os;
    os<<*chrom<<":"<<start<<"-"<<end;
    string region(os.str());
    string gene(*tmp);

    TMath::Sort(N,vals,ind);
    for (int j = 0;j < N;j++) arr[j] = vals[ind[j]];

    double median_sec = getMedianSec(arr,N);
    if (bimodal && median_sec <= 0) {
      median_sec = getMedian(arr,N);
      cout<<"Special treatment "<<gene<<" "<<region<<endl;
    }

    int n2 = 0;
    for (int i = 0;i < N;i++)
      if (arr[i] < 0.75*median_sec) break;
      else n2++;
    int n1 = N - n2;
    double median1 = getMedian(&arr[n2],n1), median2 = getMedian(arr,n2);
    double median  = getMedian(arr,N);

    double       max = 2.2*median, min = 0, range = max - min;
    if (bimodal) max = 2.2*median2;

    his->SetBins(nb,min,max);
    for (int i = 0;i < N;i++) his->Fill(arr[i]);

    double rms1  = getRMS(&arr[n2],n1,median1);
    rms1         = getRMS(&arr[n2],n1,median1,rms1);
    double rms2  = getRMS(arr,n2,median2);
    rms2         = getRMS(arr,n2,median2,rms2);
    double rms   = getRMS(arr,N,median);
    rms          = getRMS(arr,N,median,rms);
    double area1 = n1*(max - min)/nb, area2 = n2*(max - min)/nb;
    double area  = N*(max - min)/nb;

    TF1 *f1 = new TF1("fit_gaus1",gaus1,min - range,max + range,3);
    f1->SetNpx(1000);
    f1->SetParameter(0,area);
    f1->SetParameter(1,median);
    f1->SetParameter(2,rms);

    TF1 *f2 = new TF1("fit_gaus2",gaus2,min - range,max + range,6);
    f2->SetNpx(10000);
    f2->SetParameter(0,area1);
    f2->SetParameter(1,median1);
    f2->SetParameter(2,rms1);
    f2->SetParameter(3,area2);
    f2->SetParameter(4,median2);
    f2->SetParameter(5,rms2);

    TF1 *f_fit = f1; if (bimodal) f_fit = f2;
    TF1 *f_prx = (TF1*)f_fit->Clone("approx");

    his->Fit(f_fit,"qn");

    double p_prx = getKolmogorovProb1(arr,N,f_prx);
    double p_fit = getKolmogorovProb1(arr,N,f_fit);

    TF1 *func = f_fit;
    if (p_prx < p_fit) func = f_prx;

    pars[0] = func->GetParameter(0);
    pars[1] = func->GetParameter(1);
    pars[2] = func->GetParameter(2);

    if (br != NULL) br->Fill();

    delete f1;
    delete f2;
    delete f_prx;
  }

  delete his;

  tree->Write("",TObject::kOverwrite);
  file.Close();

  delete[] vals;
  delete[] arr;
  delete[] ind;
}
