#!/usr/bin/env python

import numpy as np
import argparse
from pytools import io

parser = argparse.ArgumentParser(description='Plot BAF')
parser.add_argument("root_file", help="CNVnator root file name")
parser.add_argument("-chrom", "--chromosomes", help="Comma separated chromosom list",
                    default="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22")
parser.add_argument("-bs", "--binsize", type=int,
                    help="size of bins for SNP signal in root file", default=100000)
parser.add_argument("-o", "--save_file",
                    help="save plot to file", default=None)
parser.add_argument("-t", "--title",
                    help="plot title", default=None)
parser.add_argument("-rdbs", "--rdbinsize", type=int,
                    help="size of bins for RD signal in root file if different from SNP bin size", default=100000)
parser.add_argument("-pbs", "--plotbinsize", type=int,
                    help="plot bin size", default=1000000)
parser.add_argument('-nomask', help="If not set only SNPs in P mask region will be used", action='store_true')
parser.add_argument('-useid', help="Use just SNPs that exist in given database", action='store_true')
args = parser.parse_args()

if args.save_file:
    import matplotlib as mpl

    mpl.use('Agg')
    import matplotlib.pyplot as plt
else:
    import matplotlib.pyplot as plt

bs = args.binsize
rdbs = args.rdbinsize
pbs = args.plotbinsize

iof = io.IO(args.root_file)
snpflag=0
if not args.nomask:
    snpflag|=io.FLAG_USEMASK
if args.useid:
    snpflag|=io.FLAG_USEID
print(snpflag)
chroms = iof.get_chrom_names_with_tree()
chrs = ["chr" + c if ((not c in chroms) and ("chr" + c in chroms)) else c for c in args.chromosomes.split(",")]

bafall = []
rdall = []
cstart = {}
cend = {}
Nrd = pbs / rdbs
N = pbs / bs
for c in chrs:
    snp = iof.get_signal(c, bs, "SNP maf",snpflag)[1]
    rd = iof.get_signal(c, bs, "RD")[1]
    cstart[c] = len(rdall)
    bafall += [sum(snp[n:n + N]) / len(snp[n:n + N]) for n in range(0, len(snp), N)]
    rdall += [sum(rd[n:n + N]) for n in range(0, len(snp), N)]
    cend[c] = len(rdall)

mean_rd = sum(rdall) / len(rdall)

dt = 2.0 * np.pi / len(rdall)
theta = np.arange(0, 2.0 * np.pi, dt)

fig = plt.figure()
ax = plt.subplot(111, projection='polar')

cc = (0, 0, 1)
for c in chrs:
    if cc == (0, 0, 1):
        cc = (0, 1, 0)
    else:
        cc = (0, 0, 1)
    x = []
    y = []
    for i in range(cend[c] - cstart[c]):
        if bafall[cstart[c] + i] > 0:
            x.append(np.pi / 2 - theta[cstart[c] + i])
            y.append(1. - bafall[cstart[c] + i])
    plt.polar(x, y, color=cc, linewidth=0.3)
    plt.fill_between(x, y, [1.0 for i in y], color=cc, alpha=0.8)

cc = (0.6, 0.6, 0.6)
for c in chrs:
    if cc == (0.6, 0.6, 0.6):
        cc = (0.3, 0.3, 0.3)
    else:
        cc = (0.6, 0.6, 0.6)
    x = []
    y = []
    for i in range(cend[c] - cstart[c]):
        x.append(np.pi / 2 - theta[cstart[c] + i])
        if rdall[cstart[c] + i] > (3 * mean_rd):
            y.append(1.)
        elif rdall[cstart[c] + i] < (mean_rd / 3):
            y.append(1. / 9.)
        else:
            y.append(rdall[cstart[c] + i] / (3 * mean_rd))
    plt.polar(x, y, color=cc, linewidth=0.3)
    plt.fill_between(x, [0.1 for i in y], y, color=cc, alpha=0.8)

for c in chrs:
    ax.text(np.pi / 2 - theta[(cstart[c] + cend[c]) / 2], 0.8, c, fontsize=8)

ax.set_rmax(0.9)
ax.set_rticks([])
ax.set_xticks([])
ax.grid(True)

if args.title:
    fig.suptitle(args.title, fontsize='large')
else:
    fig.suptitle(args.root_file, fontsize='large')

if args.save_file:
    plt.savefig(args.save_file, dpi=400)
else:
    plt.show()
