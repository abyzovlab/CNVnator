#include "VcfParser.hh"

VcfParser::VcfParser(string fileName) :
stdin(false),
nsamples(0),
cnames_(NULL),
clens_(NULL),
n_chr_(0),
cmax(20000),
nrec(0)
{
  // read from stdin if fileName is an empty string
  int len = fileName.length();
  if (len == 0) {
    stdin = true;
    f=bcf_open("-","r");
    if (f==NULL)
      cerr << "Vcf parser: Can't read from stdin!" << endl;
    else
      cout << "Vcf parser: reading from stdin..." << endl;
  } else {
    f=bcf_open(fileName.c_str(),"r");
    if (f==NULL) {
      cerr << "Vcf parser: Can't open file '" << fileName << "'." << endl;
      return;
    } else
      cout << "Vcf parser: reading from file '" << fileName << "'..." << endl;
  }
  clens_=new int[cmax];
  cnames_=new string[cmax];
  
  cout << "Parsing header" << endl;
  parseHeader();
  cout << "done" << endl;
  rec=bcf_init();
}

VcfParser::~VcfParser()
{
  if (cnames_) delete[] cnames_;
  if (clens_)  delete[] clens_;
  cout << "Vcf parser: finished after reading " << nrec << " records." << endl;
}

void VcfParser::parseHeader()
{
  if(f!=NULL) {
    h=bcf_hdr_read(f);
    nsamples=bcf_hdr_nsamples(h);
    n_chr_=h->n[BCF_DT_CTG];
    //cout << "Sequence names and lengths:" << endl;
    for (int i = 0; i < n_chr_; i++) {
      //cout << bcf_hdr_id2name(h, i) << " Length:" << h->id[BCF_DT_CTG][i].val->info[0] << endl;
      cnames_[i]=bcf_hdr_id2name(h, i);
      clens_[i]=h->id[BCF_DT_CTG][i].val->info[0];
    }
  }
}

int VcfParser::getChromosomeIndex() {
  int ix=0;
  while(cnames_[ix]!=c_chr) {
    ix++;
    if(ix>=n_chr_) return -1;
  }
  return ix;
}

bool VcfParser::parseRecord(bool idvar)
{
  int ngt_arr = 0;
  int ngt     = 0;
  int *gt     = NULL;
  int nad_arr = 0;
  int nad     = 0;
  int *ad     = NULL;
  bool done=false;
  char *filter=(char*)"PASS";
  if((f!=NULL)&&(h!=NULL)) {
    while(!done) {
      if(bcf_read(f,h,rec)!=0) return false;
      if (bcf_is_snp(rec) && bcf_has_filter(h,rec,filter)) {
        if(!idvar) {
          ngt = bcf_get_format_int32(h, rec, "GT", &gt, &ngt_arr);
          nad = bcf_get_format_int32(h, rec, "AD", &ad, &nad_arr);
          if(ngt==0 || nad==0) continue;
        }
        c_chr=cnames_[rec->rid];
        c_pos=rec->pos;
        c_id=rec->d.id;
        c_ref=rec->d.allele[0];
        c_qual=rec->qual;
        c_alt=rec->d.allele[1];
        c_filter="PASS";
        if(!idvar) {
          c_nref=ad[0];
          c_nalt=ad[1];
          c_gt="0/0";
          if(gt[0]==2 && gt[1]==4) c_gt="0/1";
          if(gt[0]==4 && gt[1]==4) c_gt="1/1";
          if(gt[0]==2 && gt[1]==5) c_gt="0|1";
          if(gt[0]==4 && gt[1]==5) c_gt="1|1";
          if(gt[0]==2 && gt[1]==3) c_gt="0|0";
          if(gt[0]==4 && gt[1]==3) c_gt="1|0";
        }
        nrec++;
        done=true;
      }
    }
    return true;
  } else return false;
}
