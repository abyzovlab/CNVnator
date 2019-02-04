#include "FastaParser.hh"

const int FastaParser::sbuffs=65536;

FastaParser::FastaParser(string fileName) :
buff(""),
c_name(""),
c_data(""),
n_cinb(0),
eof(false)
{
  inFile=gzopen(fileName.c_str(), "rb");
  if (inFile == NULL) {
    cerr << "Fasta parser: Can't open file '" << fileName << "'." << endl;
  }
  sbuff=new char[sbuffs];
}

FastaParser::~FastaParser()
{
  gzclose(inFile);
  delete [] sbuff;
}

bool FastaParser::nextChromosome() {
  if ((inFile==NULL) || (eof && n_cinb==0)) return false;
  while(!eof && (n_cinb<2)) {
    int n = gzread(inFile, sbuff, sbuffs-1);
    if (n>0) {
      sbuff[n]=0;
      char *pch=strchr(sbuff,'>');
      while (pch!=NULL) { n_cinb++;pch=strchr(pch+1,'>');}
      buff.append(sbuff,n);
    } else eof=true;
  }
  buff.erase(0,buff.find(">")+1);
  c_name=buff.substr(0,buff.find(" "));
  buff.erase(0,buff.find("\n")+1);
  if(!eof) {
    c_data=buff.substr(0,buff.find(">"));
    buff.erase(0,buff.find(">"));
  } else {
    c_data=buff;
    buff.clear();
  }
  n_cinb--;
  c_data.erase(std::remove(c_data.begin(), c_data.end(), '\n'),c_data.end());
  return true;
}

