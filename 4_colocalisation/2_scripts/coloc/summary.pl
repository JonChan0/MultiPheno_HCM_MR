#!/usr/bin/perl

use strict;
use warnings;

print "LOCI\tN\tH0\tH1\tH2\tH3\tH4\n";

my @t =("dbp",
	"NTproBNP",
	"bmi",
	"radialdsr");

my $sno=1;
for(my $i = 1; $i <=31; $i++){
    foreach my $t (sort(@t)){
	## bmi.hcm.LOC14.summary.txt
	open(my $F, "${t}.hcm.LOC${i}.summary.txt");

	<$F>;
	my $N = <$F>; chomp $N;
	my $H0 = <$F>; chomp $H0;
	my $H1 = <$F>; chomp $H1;
	my $H2 = <$F>; chomp $H2;
	my $H3 = <$F>; chomp $H3;
	my $H4 = <$F>; chomp $H4;

	print "LOC${i}\t$t\thcm\t$N\t$H0\t$H1\t$H2\t$H3\t$H4\n";

	close($F);
      
    }
}
