#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Plot BAF')
parser.add_argument("root_file", help="cnvnator root file name")
parser.add_argument("region", help="Chromosome or region in format chr:START-END")
parser.add_argument("-bs", "--binsize", type=int,
                    help="size of bins", default=100000)
parser.add_argument("-rdbs", "--rdbinsize", type=int,
                    help="size of bins for RD signal", default=100000)
parser.add_argument("-msc", "--minsnpc", type=int,
                    help="min number of snp-s in bin (default=10)", default=10)
parser.add_argument("-ss", "--ssize", type=int,
                    help="smalest event size in bins (default=2)", default=2)
parser.add_argument("-pmin", "--pmin", type=float,
                    help="p-threshold (0.01)", default=0.01)
parser.add_argument("-pdec", "--pdec", type=float,
                    help="p-dec", default=0.9)
parser.add_argument("-o", "--save_file",
                    help="save plot to file", default=None)
parser.add_argument("-t", "--title",
                    help="plot title", default=None)
parser.add_argument('-nomask', action='store_true')
parser.add_argument('-useid', action='store_true')
args=parser.parse_args()

pmin=-1
pmax=1e10
chr=""
bs=args.binsize
if args.region.find(":")>-1:
  ss=args.region.split(":")
  chr=ss[0]
  pmin=int(ss[1].split("-")[0])
  pmax=int(ss[1].split("-")[1])
else:
  chr=args.region


f=ROOT.TFile(args.root_file)
sn="snp_likelihood_"+chr+"_"+str(bs)
snc="snp_bafc_"+chr+"_"+str(bs)
if not args.nomask:
  sn+="_mask"
  snc+="_mask"
if args.useid:
  sn+="_id"
  snc+="_id"
fsn=f.Get("bin_"+str(bs)).Get(sn)
fsnc=f.Get("bin_"+str(bs)).Get(snc)
nsnx=fsn.GetXaxis().GetNbins()
nsny=fsn.GetYaxis().GetNbins()
m=[]
mn=[]
bafc=[]
for i in range(nsnx):
  mr=[]
  s=0
  for j in range(nsny):
    mr.append(fsn.GetBinContent(i,j))
    s+=fsn.GetBinContent(i,j)
  m.append(mr)
  mn.append(s)
  bafc.append(fsnc.GetBinContent(i))

bins=[[i] for i in range(nsnx) if bafc[i]>=args.minsnpc and not mn[i]==0.0]
lk=[m[i] for i in range(nsnx) if bafc[i]>=args.minsnpc and not mn[i]==0.0]

def pvalue(i):
  global lk,nsny
  p=0
  for k in range(nsny):
    p+=min(lk[i][k],lk[i+1][k])
  return p

def plot(lk,bins,n,iter,prefix,maxp,minp):
  global nsny,m
  mm=m[:]
  for i in range(n):
    for b in bins[i]:
      for k in range(nsny):
        mm[b][k]=lk[i][k]
  fig=plt.figure(1,figsize=(16, 9), dpi=120, facecolor='w', edgecolor='k')
  fig.suptitle("Iter: "+str(iter)+"   /   Segments: "+str(n)+"   /   Overlap interval: ("+('%.4f'%minp)+","+('%.4f'%maxp)+")", fontsize='large')
  plt.subplot(211)
  plt.ylabel("BAF")
  plt.imshow(numpy.transpose(mm),aspect='auto')
  plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
  plt.yticks([0,50.5,101,151.5,201],("1.00","0.75","0.50","0.25","0.00"))
  #plt.grid(True,color="w")
  plt.subplot(212)
  plt.xlabel("BAF")
  plt.ylabel("Likelihood")
  plt.xticks([0,0.25,0.50,0.75,1.0])
  plt.grid(True,color="b")
  for i in range(n):
    plt.plot(numpy.linspace(0,1,201),lk[i])
  plt.savefig(prefix+"_"+str(iter).zfill(4),dpi=150)
  plt.close(fig)


pv=[pvalue(i) for i in range(len(bins)-1)]  
iter=0
while len(pv)>0:
  maxp=max(pv)
  minp=max(maxp*args.pdec,args.pmin)
  #print maxp,minp,pv.index(maxp)
  if maxp<args.pmin:
    break
  i=0
  while i<(len(bins)-1):
    if pv[i]>minp:
      #print "del",i,pv[i],pv[i-1]
      nlt=[]
      ss=0
      for k in range(nsny):
        nlt.append(lk[i][k]*lk[i+1][k])
        ss=ss+lk[i][k]*lk[i+1][k]
      bins[i]+=bins[i+1]
      for k in range(nsny):
        lk[i][k]=nlt[k]/ss
      del lk[i+1]
      del bins[i+1]
      if i<(len(pv)-1):
        del pv[i+1]
      if i<(len(bins)-1):
        pv[i]=pvalue(i)
      else:
        del pv[i]
      if i>0:
        pv[i-1]=pvalue(i-1)  
      #print len(bins)
    else:
      i=i+1
  if iter%1==0:
    plot(lk,bins,len(bins),iter+1,"test0",maxp,minp)
  iter=iter+1
    
i=0
while i<len(bins):
  if len(bins[i])<args.ssize:
    del lk[i]
    del bins[i]
  else:
    i=i+1

pv=[pvalue(i) for i in range(len(bins)-1)]  
iter=0
while len(pv)>0:
  maxp=max(pv)
  minp=max(maxp*args.pdec,args.pmin)
  #print maxp,minp,pv.index(maxp)
  if maxp<args.pmin:
    break
  i=0
  while i<(len(bins)-1):
    if pv[i]>minp:
      #print "del",i,pv[i],pv[i-1]
      nlt=[]
      ss=0
      for k in range(nsny):
        nlt.append(lk[i][k]*lk[i+1][k])
        ss=ss+lk[i][k]*lk[i+1][k]
      bins[i]+=bins[i+1]
      for k in range(nsny):
        lk[i][k]=nlt[k]/ss
      del lk[i+1]
      del bins[i+1]
      if i<(len(pv)-1):
        del pv[i+1]
      if i<(len(bins)-1):
        pv[i]=pvalue(i)
      else:
        del pv[i]
      if i>0:
        pv[i-1]=pvalue(i-1)  
      #print len(bins)
    else:
      i=i+1
  if iter%1==0:
    plot(lk,bins,len(bins),iter+1,"test1",maxp,minp)
  iter=iter+1

#plt.imshow(mm,aspect='auto')
#plt.show()
#  print lk[mi]
        

exit(0)

res=args.resolution

rdbs=args.rdbinsize
rd=[]
xrd=[]

frd=f.Get("bin_"+str(rdbs)).Get("his_rd_p_"+chr+"_"+str(rdbs)+"_GC")
nrd=frd.GetSize()
for i in range(nrd):
  if ((i+1)*rdbs>pmin) and (i*rdbs<pmax):
    rd.append(frd.GetBinContent(i))
    xrd.append(i*rdbs)


fig=plt.figure(1,figsize=(12, 8), dpi=150, facecolor='w', edgecolor='k')



if args.save_file:
    plt.savefig(args.save_file,dpi=150)
    plt.close(fig)
else:
    plt.show()
