#!/usr/bin/perl

use strict;
use warnings;
use File::Basename;

my %MAP = (
    "other_allele"=>"A1",
    "effect_allele"=>"A2",
    "beta"=>"B",
    "p_value"=>"P",
    "n"=>"N",
    "N"=>"N",
    "rsid"=>"SNP",
    "n_total"=>"N");

open(my $F, "zcat $ARGV[0] |");

my $h = <$F>; chomp $h;
my @h = split "\t", $h;
my $oh = "";
my $cut = "";
for(my $i = 0; $i < scalar(@h); $i++){
    my $val = $h[$i];
    if(exists($MAP{$val})){
	$val = $MAP{$val};
	$oh.=$val."\t";
	$cut.=($i+1).",";
    }   
}
chop $oh;
chop $cut;

my $name = basename($ARGV[0]);
$name=~s/\.tsv\.gz//g; $name.=".ssimp.tsv";

open(my $OH, "> $name");
print $OH "$oh\n";
close($OH);
#print "$cut\n";


my $cmd="zcat $ARGV[0] | cut -f $cut | sed '1d' >> $name";
`$cmd`;
#print "$name\n";

close($F);
