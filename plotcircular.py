#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Plot BAF')
parser.add_argument("root_file", help="cnvnator root file name")
parser.add_argument("-chrom","--chromosomes", help="Comma separated chromosom list",default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22")
parser.add_argument("-bs", "--binsize", type=int,
                    help="size of bins", default=1000000)
parser.add_argument("-o", "--save_file",
                    help="save plot to file", default=None)
parser.add_argument("-t", "--title",
                    help="plot title", default=None)
parser.add_argument("-rdbs", "--rdbinsize", type=int,
                    help="size of bins for RD signal", default=100000)
parser.add_argument('-nomask', action='store_true')
parser.add_argument('-useid', action='store_true')

args=parser.parse_args()

chrs=args.chromosomes.split(",")
bs=args.binsize
rdbs=args.rdbinsize

f=ROOT.TFile(args.root_file)

baf={}
count={}
maxbin={}
rd={}
rdcount={}
ard=0
act=0
n=250000000/bs
for c in chrs:
  baf[c]=[0 for i in range(n)]
  count[c]=[0 for i in range(n)]
  rd[c]=[0 for i in range(n)]
  rdcount[c]=[0 for i in range(n)]
  maxbin[c]=0
  t=f.Get("vcf_"+c)
  print "Reading SNP for",c,"chromosome..."
  data=[]
  for e in t:
    if (e.nref+e.nalt)>0 and (e._gt==1 or e._gt==5 or e._gt==6) and (args.nomask or e.flag>1) and (e.flag==1 or e.flag==3 or not args.useid):
      bafv=1.0*e.nalt/(e.nref+e.nalt)
      baf[c][e.position/bs]+=(1.0-bafv) if bafv<0.5 else bafv
      count[c][e.position/bs]+=1
      if (e.position/bs)>maxbin[c]:
        maxbin[c]=e.position/bs
  frd=f.Get("bin_"+str(rdbs)).Get("his_rd_p_"+c+"_"+str(rdbs)+"_GC")
  print "Reading RD for",c,"chromosome..."
  nrd=frd.GetSize()
  for i in range(nrd):
    posb=i*rdbs/bs
    rd[c][posb]+=frd.GetBinContent(i)
    rdcount[c][posb]+=1
    ard+=frd.GetBinContent(i)
    act+=1
    if posb>maxbin[c]:
      maxbin[c]=posb

ard/=act
bafall=[]
rdall=[]
countall=0
cstart={}
cend={}
for c in chrs:
  cstart[c]=countall
  for i in range(maxbin[c]+1):
    if count[c][i]>0:
      baf[c][i]/=count[c][i]
    if rdcount[c][i]>0:
      rd[c][i]/=rdcount[c][i]
    bafall.append(baf[c][i])
    rdall.append(rd[c][i])
  countall+=maxbin[c]+1
  cend[c]=countall


dt=2.0*np.pi/countall
theta=np.arange(0,2.0*np.pi,dt)
  

ax = plt.subplot(111, projection='polar')

cc=(0,0,1)
for c in chrs:
  if cc==(0,0,1):
    cc=(0,1,0)
  else:
    cc=(0,0,1)
  x=[]
  y=[]
  for i in range(cend[c]-cstart[c]):
    if bafall[cstart[c]+i]>0:
      x.append(np.pi/2-theta[cstart[c]+i])
      y.append(bafall[cstart[c]+i])
  plt.polar(x,y,color=cc)
  plt.fill_between(x,y,[1.0 for i in y],color=cc,alpha=0.2)

cc=(0.6,0.6,0.6)
for c in chrs:
  if cc==(0.6,0.6,0.6):
    cc=(0.3,0.3,0.3)
  else:
    cc=(0.6,0.6,0.6)
  x=[]
  y=[]
  for i in range(cend[c]-cstart[c]):
    if rdall[cstart[c]+i]<(3*ard) and rdall[cstart[c]+i]>(ard/3):
      x.append(np.pi/2-theta[cstart[c]+i])
      y.append(rdall[cstart[c]+i]/(3*ard))
  plt.polar(x,y,color=cc)
  plt.fill_between(x,[0.1 for i in y],y,color=cc,alpha=0.2)

for c in chrs:
  ax.text(np.pi/2-theta[(cstart[c]+cend[c])/2], 0.8, c, fontsize=8)


ax.set_rmax(0.9)
ax.set_rticks([])
ax.set_xticks([])
ax.grid(True)


if args.save_file:
    plt.savefig(args.save_file,dpi=400)
else:
    plt.show()
