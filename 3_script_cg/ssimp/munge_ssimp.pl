#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $inf = $ARGV[2];

my %HASH = ();

open(my $F1, $ARGV[1]);

while(my $l = <$F1>){
    chomp $l;
    my ($chr, $pos0, $pos, $snp) = split "\t", $l;
    $chr=~s/chr//g;
    my %subhash = ();    
    $subhash{"chr"}=$chr;
    $subhash{"pos"}=$pos;

    push @{$HASH{$snp}}, \%subhash;
}

close($F1);

#print Dumper(%HASH);

## fields required for hyprcoloc
# chromosome
# base_pair_location
# effect_allele
# other_allele
# beta
# standard_error

open(my $F2, $ARGV[0]);

my $h = <$F2>; chomp $h;

print "chromosome\tbase_pair_location\teffect_allele\tother_allele\tbeta\tstandard_error\tp_value\trsid\tsource\tinfo\tbp_37\n";

while(my $l = <$F2>){
     chomp $l;
     my @d = split "\t", $l;
     my $snp = $d[4];

     if(exists($HASH{$snp})){
	 if(scalar(@{$HASH{$snp}}) > 1){
	     print STDERR "Mutiple pos for SNP $snp\n";
	 }
	 
	 my $chr=$d[0];
	 my $bp_37 = $d[1];
	 my $subhash = $HASH{$snp}[0];
	 
	 if($chr != $$subhash{"chr"}){    print STDERR "chromosome mismatch for SNP $snp\n";next;}
	 
	 my $pos = $$subhash{"pos"};
#	 print "$l\t$pos\n";

	 my $A1 = $d[5]; my $A2 = $d[6];
	 my $beta = $d[14];
	 ## calculate se
	 my $N = $d[13];

	 if($N==0){next;}
	 
	 my $se = 1/(sqrt($N));
	 my $z = $d[2]; my $p = $d[12];

	 my $source = $d[3]; my $info = $d[8];

	 # 1       chromosome
	 # 2       base_pair_location
	 # 3       effect_allele
	 # 4       other_allele
	 # 5       beta
	 # 6       standard_error
	 # 8       p_value
	 # 9       rsid


	 if($info>=$inf && $source eq "SSIMP"){	 
	     print "$chr\t$pos\t$A2\t$A1\t$beta\t$se\t$p\t$snp\t$source\t$info\t$bp_37\n";
	 }
#	 print "$z\t$p\n";
#	 last;
     }
}

close($F2);
