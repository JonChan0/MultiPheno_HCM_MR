use strict;
use warnings;
use Data::Dumper;

my $chrin=$ARGV[0]; my $st=$ARGV[1]; my $en=$ARGV[2]; ## coords
my $gwas_file=$ARGV[3];
my $out_prefix=$ARGV[4];

my %HASH = ();
my %SNPHASH = ();
my @PHENO = ();

open(my $F, $gwas_file);
open(my $LOG, "> $out_prefix.log");

while(my $l = <$F>){
    chomp $l;
    if($l=~/^#/){next;}

    my @d = split "\t", $l;
    my $pheno = $d[1];
    my $file = $d[2];
    
    print $LOG "$pheno $file\n";
    print STDERR "$pheno $file\n";

    push @PHENO, $pheno;

    open(my $F1, "zcat DATASETS/$file |");

    my $h = <$F1>; chomp $h;
    my @h = split "\t", $h; my %h = ();
    for(my $i = 0; $i < scalar(@h); $i++){
	$h{$h[$i]} = $i;
    }

    my $count = 0;
    while(my $l2 = <$F1>){
	chomp $l2;
	my @d = split "\t", $l2;
	my $chr = $d[$h{"chromosome"}];
	my $bp = $d[$h{"base_pair_location"}];

	if($chr eq "NA") {next;}

	my $a1 = $d[$h{"effect_allele"}];
	my $a2 = $d[$h{"other_allele"}];

	my $beta = $d[$h{"beta"}];
	my $se = $d[$h{"standard_error"}];

	if($beta eq "NA" || $se eq "NA"){
	    next;
	}
	
	# print "$snp\t$chr\t$bp\t$a1\t$a2\t$beta\t$se\n";

	if($chr eq "X"){next;}
	
	if($chr==$chrin && $bp >= $st && $bp <= $en){
	    my $snptest1=$chr.":".$bp."_".$a1."_".$a2;
	    my $snptest2=$chr.":".$bp."_".$a2."_".$a1;

	    my $snpused=$snptest1;
	
	    if(exists($SNPHASH{$bp}{$snptest2})){
		$snpused=$snptest2;
	    }
	    elsif(!exists($SNPHASH{$bp}{$snptest1})){
		$SNPHASH{$bp}{$snptest1}++;
	    }

	    $HASH{$snpused}{$pheno}{"a1"} = $a1;
	    $HASH{$snpused}{$pheno}{"beta"} = $beta;
	    $HASH{$snpused}{$pheno}{"se"} = $se;
	    $count++;
	}
    }
    
    close($F1);

#    last;

    print $LOG "$pheno: $count\n"; ## log the number of variants per phenotype

#    last;
}

close($F);

#print Dumper(%HASH);

#exit;

open(my $OBETA, "> $out_prefix.beta.txt");
open(my $OSE, "> $out_prefix.se.txt");

my $hl = "snp\t";
foreach my $pheno (@PHENO){
    $hl=$hl."$pheno\t";
}
chop $hl;
print $OBETA "$hl\n";
print $OSE "$hl\n";

my $totlen=0;
foreach my $bp (sort{$a<=>$b}(keys(%SNPHASH))){
    foreach my $snptestid (keys(%{$SNPHASH{$bp}})){
	my $betal="$snptestid\t";
	my $sel="$snptestid\t";

	my $alignedallele="NA";
	
	foreach my $pheno (@PHENO){
	    if(!exists($HASH{$snptestid}{$pheno})){ ## if the SNP is not present for one or more phenotypes.
		my $failedlist = "";
		my $failct=0;

		foreach my $pheno (@PHENO){
		    if(! exists($HASH{$snptestid}{$pheno})){
			$failedlist=$failedlist.";".$pheno;
			$failct++;
		    }
		}
		print $LOG "$snptestid failed: $failct $failedlist\n";
		goto END;
	    }

	    my $beta = $HASH{$snptestid}{$pheno}{"beta"};
	    my $se = $HASH{$snptestid}{$pheno}{"se"};

	    my $a1 = $HASH{$snptestid}{$pheno}{"a1"};

	    if($alignedallele eq "NA"){$alignedallele=$a1;}

	    if($alignedallele ne $a1){$beta=$beta*-1;}
	    
	    $betal=$betal.$beta."\t";
	    $sel=$sel.$se."\t";
	}

	chop $betal;
	print $OBETA "$betal\n";

	chop $sel;
	print $OSE "$sel\n";
	
	$totlen++;
      END:
    }
}

print $LOG "Total: $totlen\n"; ## total number of SNPs.

close($OBETA);
close($OSE);
close($LOG);
