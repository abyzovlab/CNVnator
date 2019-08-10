"""CNVnator python tools::
    
    ROOT IO class
"""

from __future__ import print_function
import ROOT
import sys

FLAG_SEX = 0x0001
FLAG_GC_CORR = 0x0010
FLAG_AT_CORR = 0x0020
FLAG_USEMASK = 0x0100
FLAG_USEID = 0x0200
FLAG_USEHAP = 0x0400


class IO:
    signals = {
        "RD": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD unique": "his_rd_u_%(chr)s_%(bin_size)d%(rd_flag)s",
        "RD raw": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_raw",
        "RD l1": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l1",
        "RD l2": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l2",
        "RD l3": "his_rd_p_%(chr)s_%(bin_size)d%(rd_flag)s_l3",
        "RD partition": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s",
        "RD call": "his_rd_p_%(chr)s_%(bin_size)d_partition%(rd_flag)s_merge",
        "GC": "%(chr)s_gc_%(bin_size)",
        "SNP count": "snp_bafc_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP baf": "snp_baf_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP maf": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP likelihood": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i1": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i2": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i3": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP i4": "snp_i4_bafc_%(chr)s_%(bin_size)d%(snp_flag)s",
        "SNP likelihood partition": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP maf partition": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i1 partition": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i2 partition": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i3 partition": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP i4 partition": "snp_i4_%(chr)s_%(bin_size)d%(snp_flag)s_partition",
        "SNP likelihood call": "snp_likelihood_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP maf call": "snp_maf_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i1 call": "snp_i1_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i2 call": "snp_i2_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i3 call": "snp_i3_%(chr)s_%(bin_size)d%(snp_flag)s_call",
        "SNP i4 call": "snp_i4_%(chr)s_%(bin_size)d%(snp_flag)s_call",
    }

    def __init__(self, root_file):
        """Class constructor
        Opens root file"""
        try:
            self.file = ROOT.TFile(root_file)
        except IOError as ioerror:
            print("File not found!")

    def __del__(self):
        """Class destructor
        Closes root file"""
        self.file.Close()

    def sufix_rd_flag(self, flags):
        """ Converts binary flags into sufix used in RD signal names """
        s = ""
        if flags & FLAG_AT_CORR:
            s += "_AT"
        if flags & FLAG_GC_CORR:
            s += "_GC"
        return s

    def sufix_snp_flag(self, flags):
        """ Converts binary flags into sufix used in SNP signal names """
        s = ""
        if flags & FLAG_USEMASK:
            s += "_mask"
        if flags & FLAG_USEID:
            s += "_id"
        if flags & FLAG_USEHAP:
            s += "_hap"
        return s

    def sufix_flag(self, flags):
        """ Converts binary flags into sufix used in distribution names """
        s = ""
        if flags & FLAG_SEX:
            s += "_sex"
        return s

    def treeName(self, chr, signal):
        """Returns TTree name for chromosome"""
        if signal == "RD":
            return chr
        else:
            return "snp_%s" % chr

    def signal_name(self, chr, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """Returns TH1 or TH2 name for signal"""
        return self.signals[signal] % {"chr": chr, "bin_size": bin_size, "rd_flag": self.sufix_rd_flag(flags),
                                       "snp_flag": self.sufix_snp_flag(flags), "flag": self.sufix_flag(flags)}

    def get_chrom_names_with_tree(self,snp=False):
        iter = self.file.GetListOfKeys();
        forest=[]
        for i in iter:
            if i.GetClassName()=="TTree":
                name=i.GetName()
                if snp and name.find("snp_") == 0:
                    forest.append(name[4:])
                if (not snp) and name.find("snp_") != 0:
                    forest.append(name)
        return forest


    def get_tree(self, chr, signal):
        """ToDo - read tree and return arrays"""
        return True

    def get_signal(self, chr, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """Returns 1D histogram: (Xbin_centers, Values) """
        his = self.file.Get("bin_" + str(bin_size)).Get(self.signal_name(chr, bin_size, signal, flags))
        n = his.GetSize()
        return [his.GetBinCenter(i) for i in range(n)], [his.GetBinContent(i) for i in range(n)]

    def get_signal_2d(self, chr, bin_size, signal, flags=FLAG_USEMASK | FLAG_GC_CORR):
        """ToDo - Returns 2D histogram: (Xbin_centers, Values) """
        return True


if __name__ == '__main__':
    print("main")
    io = IO(sys.argv[1])
    print(io.get_chrom_names_with_tree(snp=True))
