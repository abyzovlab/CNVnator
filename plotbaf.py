#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy
import ROOT
import argparse

parser = argparse.ArgumentParser(description='Plot BAF')
parser.add_argument("root_file", help="cnvnator root file name")
parser.add_argument("region", help="Chromosome or region in format chr:START-END")
parser.add_argument("-bs", "--binsize", type=int,
                    help="size of bins", default=100000)
parser.add_argument("-res", "--resolution", type=int,
                    help="size of bins", default=100)
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
if args.region.find(":")>-1:
  ss=args.region.split(":")
  chr=ss[0]
  pmin=int(ss[1].split("-")[0])
  pmax=int(ss[1].split("-")[1])
else:
  chr=args.region


f=ROOT.TFile(args.root_file)
t=f.Get("snp_"+chr)
data=[]
for e in t:
  if e.position>pmin and e.position<pmax:
    data.append([e.position,e._gt,e.nref,e.nalt,e.flag])

res=args.resolution
bs=args.binsize

def f(k,m,p):
  if(k==m):
    return 1.0*p**k*(1-p)**m
  else:
    return 1.0*p**k*(1-p)**m+1.0*p**m*(1-p)**k
    
def lh(samples):
  if len(samples)==0:
    return numpy.zeros(res+1)
  x=numpy.arange(0,1.+1.0/res,1.0/res)
  y=numpy.ones(res+1)    
  for (a,b) in samples:
    s=0.
    for pi in range(res+1):
      y[pi]=y[pi]*f(a,b,x[pi])
      s=s+y[pi]
    y=y/s
  return y
  
g1=[]
g2=[]
baf=[]
pos=[]
for i in data:
  if (i[3]+i[2])>0 and (i[1]==1 or i[1]==5 or i[1]==6) and (args.nomask or i[4]>1) and (i[4]==1 or i[4]==3 or not args.useid):
    pos.append(i[0])
    nr=i[2]
    na=i[3]
    if na<nr:
      na=na+1
    elif nr<na:
      nr=nr+1
    g1.append(nr)
    g2.append(na)    
    baf.append(1.0*na/(nr+na))

def maf(x):
  if x>0.5:
    return 1-x
  else:
    return x

n=(pos[-1]-pos[0])/bs+1

binp=[pos[0]+i*bs for i in range(n)]
bg=[[] for i in range(n)]
bbaf=[[] for i in range(n)]
bcount=[0 for i in range(n)]

for i in range(len(baf)):
  bix=(pos[i]-pos[0])/bs
  bg[bix].append((g1[i],g2[i]))
  bbaf[bix].append(baf[i])
  bcount[bix]=bcount[bix]+1


m=numpy.vstack([lh(bg[i]) for i in range(n)])
av=numpy.array([0.5 if bcount[i]==0 else numpy.mean(map(maf,bbaf[i])) for i in range(n)])
sd=numpy.array([0 if bcount[i]==0 else numpy.std(map(maf,bbaf[i])) for i in range(n)])

ng=[]
ng2=[]
for i in range(n):
  m1=0
  m2=0
  i1=0
  i2=0
  for j in range(res+1):
    if j>res/2:
      if m[i][j]>m1:
        m1=m[i][j]
        i1=j
    elif j==res/2:
      if m[i][j]>m1:
        m1=m[i][j]
        i1=j
      if m[i][j]>m2:
        m2=m[i][j]
        i2=j
    else:
      if m[i][j]>m2:
        m2=m[i][j]
        i2=j
  ng.append(1.0*(i1-i2)/res)
  ng2.append(1.0*m[i][res/2]/m1)  

fig=plt.figure(1,figsize=(12, 8), dpi=150, facecolor='w', edgecolor='k')
plt.subplot(411)
plt.ylabel("BAF")
plt.scatter(pos,baf,s=2,c="b", alpha=0.5,marker='.')
plt.axis([pos[0], pos[-1], 0, 1])
plt.yticks([0,0.25,0.5,0.75,1.0],("0","0.25","0.50","0.75","1.00"))
plt.grid(True)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)

plt.subplot(412)
plt.ylabel("Likelihood")
plt.imshow(numpy.transpose(m),aspect='auto')
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
plt.yticks([0,res/4,res/2-1,3*res/4-1,res+1],("1.00","0.75","0.50","0.25","0.00"))
plt.grid(True,color="w")

plt.subplot(413)
plt.ylabel("W, C/M")

plt.plot(binp,ng,'r-',binp,ng2,'b.', markersize=5)
plt.axis([binp[0], binp[-1], -0.05, 1.05])
plt.yticks([0,0.25,0.5,0.75,1.0],("0","0.25","0.50","0.75","1.00"))
plt.grid(True)
plt.tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)

plt.subplot(414)
plt.ylabel("MAF")
plt.errorbar(binp, av, yerr=sd, fmt='o',marker='o', mfc='red', mec='green', ms=0.1, mew=0.0)
plt.plot(binp, av,alpha=0.5,marker='.')
plt.axis([binp[0], binp[-1], 0, 0.6])
plt.yticks([0,0.25,0.5,0.75,1.0],("0","0.25","0.50","0.75","1.00"))
plt.grid(True)
if args.title:
  fig.suptitle(args.title, fontsize='large')
else:
  fig.suptitle(args.root_file+" - "+chr+":"+str(pos[0])+"-"+str(pos[-1]), fontsize='large')

if args.save_file:
    plt.savefig(args.save_file,dpi=150)
    plt.close(fig)
else:
    plt.show()
