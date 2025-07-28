#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $label = $ARGV[2];

my %LOCI = ();

open(my $F0, "../../hcm_nonMTAG_31loci_b38.txt");

while(my $l = <$F0>){
    chomp $l;
    my ($sno, $snp, $chr, $st, $en, $loci) = split "\t", $l;
    my %subhash = ();
    $subhash{"chr"}=$chr;
    $subhash{"st"}=$st;
    $subhash{"en"}=$en;
    $subhash{"loci"}=$loci;

    my $file = $label.".".$loci.".txt";
    open(my $O, "> $file");
    print $O "beta\tvarbeta\tsnp\tposition\ttype\tN\tMAF\n";
    $subhash{"file"}=$O;

    push @{$LOCI{$chr}}, \%subhash;

    #last;
}

close($F0);

#print Dumper(%LOCI);

#exit;

my $type = $ARGV[1];
open(my $F1, "zcat $ARGV[0] |");

my $h = <$F1>; chomp $h;
my @h = split "\t", $h;
my %h = ();
for(my $i = 0; $i < scalar(@h); $i++){
    $h{$h[$i]}=$i;
}

#print Dumper(%h);

while(my $l = <$F1>){
    chomp $l;
    my @d = split "\t", $l;

    my $beta = $d[$h{"beta"}];
    my $se = $d[$h{"standard_error"}];

    if($se eq "NA"){next;}
    
    my $varbeta = $se**2;

    my $snp = $d[$h{"rsid"}];
    my $position = $d[$h{"base_pair_location"}]; ## build 38

    my $n = $d[$h{"n"}];
    my $eaf = $d[$h{"effect_allele_frequency"}];
    my $maf = $eaf;
    if($maf > 0.5){$maf=1-$maf;}

    my $chr = $d[$h{"chromosome"}];
    
    foreach my $subhash (@{$LOCI{$chr}}){
	my $st = $$subhash{"st"};
	my $en = $$subhash{"en"};

	if($position >= $st && $position <= $en){
	    my $O = $$subhash{"file"};
	    print $O "$beta\t$varbeta\t$snp\t$position\t$type\t$n\t$maf\n";	    
	}
    }   
#    last;
}

close($F1);

foreach my $chr (keys(%LOCI)){
    foreach my $subhash (@{$LOCI{$chr}}){
	my $O = $$subhash{"file"};
	close($O);
    }
}

