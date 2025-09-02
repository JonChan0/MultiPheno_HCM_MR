#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $DATASETS="../DATASETS/MR_DATA_nonMTAG/";

my %HASH = ();

my $file = $ARGV[0];
my $prefix = $file; $prefix=~s/\.tsv\.gz//g;

open(my $F0, "zcat ${DATASETS}${file} |");
open(my $O, "> ${DATASETS}${prefix}.with_ssimp.tsv");

print STDERR "GWAS data...\n";

## chromosome      base_pair_location      effect_allele   other_allele    beta   standard_error   p_value rsid    source  info
my $h = <$F0>; chomp $h;
my @h = split "\t", $h;
my %h = ();
for(my $i = 0; $i < scalar(@h); $i++){
    $h{$h[$i]}=$i;
}

print $O "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value\trsid\tsource\tinfo\tbp_37\n";

while(my $l = <$F0>){
    chomp $l;
    my @d = split "\t", $l;

    print $O $d[$h{"chromosome"}]."\t".$d[$h{"base_pair_location"}]."\t".$d[$h{"effect_allele"}]."\t".$d[$h{"other_allele"}]."\t";
    print $O $d[$h{"beta"}]."\t".$d[$h{"standard_error"}]."\t".$d[$h{"p_value"}]."\t".$d[$h{"rsid"}]."\t"."GWAS\t1\tNA";
    print $O "\n";
}

close($F0);

print STDERR "SSIMP data...\n";

opendir(my $D, "chunks_mr_nonMTAG_output");

my @files = grep(/$prefix/,readdir($D));

foreach my $file (@files){
    open(my $F1, "chunks_mr_nonMTAG_output/".$file);

    my $h = <$F1>; chomp $h;

    while(my $l = <$F1>){
	chomp $l;
	print $O "$l\n";
    }
    
    close($F1);
}

close($D);
close($O);
