#!/usr/bin/perl

use warnings;
use strict;
use Getopt::Long;

my $usage = "\tcnvnator2VCF.pl [-prefix prefix] file.calls [genome_dir]\n";

my ($file,$dir,$prefix);
die("not enough argyments. $usage\n") unless ( @ARGV );
GetOptions( 'p|prefix:s' => \$prefix);
$file = shift @ARGV;
if( @ARGV ) {
    $dir = shift @ARGV;
} else {
    $dir = "./";
}

if (! defined $file || ! -f $file ) {
    print STDERR $usage,"\n";
    exit;
}

open(FILE,$file) or die "Can't open file ",$file,".\n";
print STDERR "Reading calls ...\n";
my ($pop_id) = split(/\./,$file);
print '##fileformat=VCFv4.1',"\n";
print '##fileDate='.`date '+%Y%m%d'`;
print '##reference=1000GenomesPhase3_decoy-GRCh37',"\n";
print '##source=CNVnator',"\n";
print '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',"\n";
print '##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">',"\n";
print '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">',"\n";
print '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',"\n";
print '##INFO=<ID=natorRD,Number=1,Type=Float,Description="Normalized RD">',"\n";
print '##INFO=<ID=natorP1,Number=1,Type=Float,Description="e-val by t-test">',"\n";
print '##INFO=<ID=natorP2,Number=1,Type=Float,Description="e-val by Gaussian tail">',"\n";
print '##INFO=<ID=natorP3,Number=1,Type=Float,Description="e-val by t-test (middle)">',"\n";
print '##INFO=<ID=natorP4,Number=1,Type=Float,Description="e-val by Gaussian tail (middle)">',"\n";
print '##INFO=<ID=natorQ0,Number=1,Type=Float,Description="Fraction of reads with 0 mapping quality">',"\n";
print '##INFO=<ID=natorPE,Number=1,Type=Integer,Description="Number of paired-ends support the event">',"\n";
print '##INFO=<ID=SAMPLES,Number=.,Type=String,Description="Sample genotyped to have the variant">',"\n";
print '##ALT=<ID=DEL,Description="Deletion">',"\n";
print '##ALT=<ID=DUP,Description="Duplication">',"\n";
print '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',"\n";
print '##FORMAT=<ID=CN,Number=1,Type=Integer,Description="Copy number genotype for imprecise events">',"\n";
print '##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Number of paired-ends that support the event">',"\n";
print "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$pop_id\n";
my ($prev_chrom,$chrom_seq,$count) = ("","",0);
while (my $line = <FILE>) {
    my ($type,$coor,$len,$rd,$p1,$p2,$p3,$p4,$q0,$pe) = split(/\s+/,$line);
    my ($chrom,$start,$end) = split(/[\:\-]/,$coor);
    my $isDel = ($type eq "deletion");
    my $isDup = ($type eq "duplication");
    if ($isDup) {
    } elsif ($isDel) {
    } else {
	print STDERR "Skipping unrecognized event type '",$type,"'.\n";
	next;
    }
    if ($chrom ne $prev_chrom) {
	$chrom_seq  = parseChrom($chrom,$dir);
	$prev_chrom = $chrom;
    }
    $count++;
    my $id = "";
    if ( defined $prefix ) { $id = $prefix."_"; }
    $id .= "CNVnator_";
    if    ($isDel) { $id .= "del_"; }
    elsif ($isDup) { $id .= "dup_"; }
    $id .= $count;
    print $chrom,"\t",$start,"\t",$id,"\t";
    if ($start < length($chrom_seq)) {
	print substr($chrom_seq,$start - 1,1),"\t";
    } else {
	print "N\t";
    }
    if    ($isDel) { print "<DEL>"; }
    elsif ($isDup) { print "<DUP>"; }
    print "\t.\tPASS\t";
    my $INFO = "END=".$end;
    if    ($isDel) {
	$INFO .= ";SVTYPE=DEL";
	$INFO .= ";SVLEN=-".int($len);
    } elsif ($isDup) {
	$INFO .= ";SVTYPE=DUP";
	$INFO .= ";SVLEN=".int($len);
    }
    $INFO   .= ";IMPRECISE";
    if (defined($rd) && ($rd ne "")) { $INFO .= ";natorRD=".$rd; }
    if (defined($p1) && ($p1 ne "")) { $INFO .= ";natorP1=".$p1; }
    if (defined($p2) && ($p2 ne "")) { $INFO .= ";natorP2=".$p2; }
    if (defined($p3) && ($p3 ne "")) { $INFO .= ";natorP3=".$p3; }
    if (defined($p4) && ($p4 ne "")) { $INFO .= ";natorP4=".$p4; }
    if (defined($q0) && ($q0 ne "")) { $INFO .= ";natorQ0=".$q0; }
    if (defined($pe) && ($pe ne "")) { $INFO .= ";natorPE=".$pe; }
    print $INFO;

    my $GT="GT";

    if(defined($rd) && ($rd ne "")) {	
	$GT.=":CN";
	if(defined($pe) && ($pe ne "")) {
	    $GT.=":PE";
	}
	$GT.="\t";

	if ($isDel && $rd < 0.25) {
	    $GT .= "1/1:0";   
	} elsif ($isDel && $rd >= 0.25) {
	    $GT .= "0/1:1";
	} elsif ($isDup && $rd <= 1.75) {
	    $GT .= "0/1:2";
	} elsif ($isDup && $rd > 1.75 && $rd <= 2.25) {
	    $GT .= "1/1:2";
	} elsif ($isDup && $rd > 2.25) {
	    my $cn = sprintf("%.0f",$rd);
	    $GT.="./1:$cn"; # w/o other data, we can't really say if this is
	                    # a hom dup, or het dup with higher copy number.
	} else {
	    $GT = "GT\t./.";
	}

	if(defined($pe) && ($pe ne "")) {
	    $GT.=":$pe";
	}
    } else {
	$GT.="\t./.";
    }
    print "\t$GT\n";
}
close(FILE);

exit;

sub parseChrom
{
    my ($chrom,$dir) = @_;
    my $file = $dir."/".$chrom.".fa";
    if (!open(SEQFILE,"$file")) {
        print STDERR "Can't parse sequence for chromosome ",$chrom,".\n";
        return "";
    }
    my ($header,$seq) = ("","");
    while (my $line = <SEQFILE>) {
        chomp($line);
        if (length($line) <= 0) { next; }
        if (substr($line,0,1) eq ">") {
            $header = substr($line,1);
            $seq = "";
        } else { $seq .= $line; }
    }
    close(SEQFILE);
    return $seq;
}
