#!/usr/bin/python
import ROOT
import argparse

parser = argparse.ArgumentParser(description='List content of CNVnator file')
parser.add_argument("root_file", help="cnvnator root file name")
args=parser.parse_args()

f=ROOT.TFile(args.root_file)
k=f.GetListOfKeys()
cont=[i.GetName() for i in k if i.IsFolder()]
rd=[i for i in cont if i.find("bin_")!=0 and i.find("vcf_")!=0]
vcf=[i for i in cont if i.find("vcf_")==0]
bin=[i for i in cont if i.find("bin_")==0]

print "File",args.root_file,"contains following data:"
print
ix=0
if len(rd)>0:
  ix=ix+1
  print str(ix)+") Read mapping for chromosomes:"
  print ", ".join(rd)
  print
if len(vcf)>0:
  ix=ix+1
  print str(ix)+") SNPs for chromosomes:"
  print ", ".join(map(lambda x:x[4:],vcf))
  print
if len(bin)>0:
  ix=ix+1
  print str(ix)+") Folders with histograms for bin sizes:"
  print ", ".join(map(lambda x:x[4:],bin))
  print

if ix==0:
  print "No CNVnator data!"
