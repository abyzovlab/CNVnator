#include "HisMaker.hh"
#include "Genotyper.hh"

Genotyper::Genotyper(HisMaker *maker,
		     TString file,
		     int bin) : _maker(maker),
				_file(file),
				_bin(bin),
				_hisSignal(NULL),
				_hisDistr(NULL),_hisDistrAll(NULL),
				_hisDistr1000(NULL),_hisDistr1000All(NULL)
						     
{}

void Genotyper::printGenotype(string chr,int start,int end,
			      bool useATcorr,bool useGCcorr)
{
  if (!_maker) return;

  if (start == 0) start++;
  if (start < 0 || end < 0 || end < start) {
    cerr<<"Bad coordinates: "<<chr<<" "<<start<<" "<<end<<endl;
    return;
  }

  string alt_chr              = Genome::canonicalChromName(chr);
  if (alt_chr == chr) alt_chr = Genome::extendedChromName(chr);
  TString nameSignal        = _maker->getSignalName(chr,  _bin,
						    useATcorr,useGCcorr);
  TString alt_nameSignal    = _maker->getSignalName(alt_chr,_bin,
						    useATcorr,useGCcorr);
  TString nameDistr         = _maker->getDistrName(chr,   _bin,
						   useATcorr,useGCcorr);
  TString alt_nameDistr     = _maker->getDistrName(alt_chr, _bin,
						   useATcorr,useGCcorr);
  TString nameDistrAll      = _maker->getDistrName(Genome::CHRALL,_bin,
						   useATcorr,useGCcorr);
  TString nameDistr1000     = _maker->getDistrName(chr,   1000,
						   useATcorr,useGCcorr);
  TString alt_nameDistr1000 = _maker->getDistrName(alt_chr, 1000,
						   useATcorr,useGCcorr);
  TString nameDistr1000All  = _maker->getDistrName(Genome::CHRALL,1000,
						   useATcorr,useGCcorr);
  
  TString dirName     = _maker->getDirName(_bin);
  TString dirName1000 = _maker->getDirName(1000);

  if (!_hisSignal || nameSignal != _hisSignal->GetName())
    _hisSignal = _maker->getHistogram(_file,dirName,nameSignal,alt_nameSignal);
  
  if (!_hisDistr || nameDistr != _hisDistr->GetName()) {
    _hisDistr  = _maker->getHistogram(_file,dirName,nameDistr,alt_nameDistr);
    if (_hisDistr)
      _maker->getMeanSigma(_hisDistr,_mean,_sigma);
  }
  
  if (!_hisDistrAll || nameDistrAll != _hisDistrAll->GetName()) {
    _hisDistrAll = _maker->getHistogram(_file,dirName,nameDistrAll);
    if (_hisDistrAll)
      _maker->getMeanSigma(_hisDistrAll,_meanAll,_sigmaAll);
  }

  if (!_hisDistr1000 || nameDistr1000 != _hisDistr1000->GetName()) {
    _hisDistr1000 = _maker->getHistogram(_file,dirName1000,
					 nameDistr1000,alt_nameDistr1000);
    if (_hisDistr1000)
      _maker->getMeanSigma(_hisDistr1000,_mean1000,_sigma1000);
  }
  
  if (!_hisDistr1000All || nameDistr1000All != _hisDistr1000All->GetName()) {
    _hisDistr1000All = _maker->getHistogram(_file,dirName1000,
					    nameDistr1000All);
    if (_hisDistr1000All)
      _maker->getMeanSigma(_hisDistr1000All,_mean1000All,_sigma1000All);
  }

  double scale = 2;
  if (Genome::isSexChrom(chr))
    if (_meanAll > 0 && _mean/_meanAll < 0.66) {
      cout<<"Assuming male individual!"<<endl;
      scale = 1;
    }

  double gen1 = -1,gen2 = -1;
  double av1 = (end - start + 1)*_mean/_bin;
  double av2 = (end - start + 1)*_mean1000/1000;
  if (_hisSignal && _hisDistr && _mean > 1)
    gen1 = scale*getReadCount(start,end)/av1;
  if (_hisDistr1000 && _mean1000 > 1)
    gen2 = scale*getReadCount(start,end)/av2;

  if (_hisSignal && _hisDistr)
    cout<<"Genotype "<<chr<<":"<<start<<"-"<<end<<" "<<_file<<" "
	<<gen1<<" "<<gen2<<endl;
  else
    cerr<<"Can't find all necessary histograms for chromosome '"
	<<chr<<"' in file "<<_file<<"."<<endl;
}

double Genotyper::getReadCount(int start,int end)
{
  double signal = 0, bin_over = 1./_bin;
  int bin_start = int(start*bin_over) + 1;
  int bin_end   = int(end*bin_over) + 1;
  if (bin_start == bin_end) {
    double fr = (end - start + 1)*bin_over;
    signal = _hisSignal->GetBinContent(bin_start)*fr;
  } else {
    double fr_start = (bin_start*_bin - start)*bin_over;
    double fr_end   = (end - (bin_end - 1)*_bin)*bin_over;
    signal += _hisSignal->GetBinContent(bin_start)*fr_start;
    signal += _hisSignal->GetBinContent(bin_end)*fr_end;
    for (int i = bin_start + 1;i < bin_end;i++)
      signal += _hisSignal->GetBinContent(i);
  }
  return signal;
}
