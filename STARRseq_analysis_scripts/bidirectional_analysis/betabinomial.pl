#!/usr/bin/perl
use strict;
use warnings;
use Sort::Naturally;
use Math::CDF qw(:all);
use Math::Random qw(:all);
use PDL;
use PDL::GSLSF::GAMMA;
use PDL::Stats::Basic;

die "usage: $0 <dACRfragments.bed>\n" unless @ARGV == 1;

print "#dACRchr\tdACRstart\tdACRend\tdACRactivity\tFRAGchr\tFRAGstart\tFRAGend\tRNA_F\tRNA_R\tINPUT_F\tINPUT_R\tEmpirical_bias\tSim_bias\tBinomialTest\tProb\tSim_Prob\n";


## read files
print STDERR "iterating over sites...\n";
my $counts = 0;
open F, $ARGV[0] or die;

# iterate
while(<F>){
	$counts++;
	chomp;

        # print update
        if(($counts % 1000000) == 0){
                print STDERR "##### scanned $counts sites #####\n";
        }	

	# arrays
	my @col = split("\t",$_);

	# skip if read depths are too low
	if(($col[7]+$col[8]+$col[9]+$col[10]) == 0){
		next;
	}

	# estimate beta-binomial
	else{	

		my $a1 = $col[7];
		my $a2 = $col[8];
		my $at = $a1+$a2;
		my $b1 = $col[9];
		my $b2 = $col[10];
		if($a1 > 0 && $b1 == 0){
			$b1 = $b1 + 1;
		}
		if($a2 > 0 && $b2 == 0){
			$b2 = $b2 + 1;
		}
		my $bt = $b1+$b2;
		my $pa = ($a1 + $b1)/($a1 + $a2 + $b1 + $b2);
		my $p1 = $a1/($at);
		my $p2 = $b1/($bt);
		my $pb = $b2/($bt);
		my $num = $p1 - $p2;
		if($pa == 0){
			$pa = 0.1
		}
		my $den1 = sqrt(($pa*(1-$pa))*((1/$at)+(1/$bt)));
		my $prob = 1;
		if($den1 == 0){
			$prob = 1;
		}
		else{
			my $ts = $num/$den1;
			my $score = - abs($ts);
			$prob = 2*(pnorm($score));
		}
		my @bb = betabinom(10000, $a1, $a2, $b1, $b2);
		my $binomtest;
		if($p2 > 0){
			my $x = pdl($a1);
			my $n = pdl($at);
			my $p = pdl($p2);
			#print STDERR "binomial: x=$a1 n=$at p=$p2\n";
			$binomtest = binomial_test($x, $n, $p);
		}
		else{
			my $x = pdl($a2);
			my $n = pdl($at);
			my $p = pdl($pb);
			#print STDERR "binomial: x=$a2 n=$at p=$pb\n";
			$binomtest = binomial_test($x, $n, $p);
		}
		my $dif = abs(($a1/($a1+$a2))-($b1/($b1+$b2)));
		$dif = sprintf "%.3f", $dif;
		$bb[1] = sprintf "%.3f", $bb[1];
		print "$_\t$dif\t$bb[1]\t$binomtest\t$prob\t$bb[0]\n";
		next;
	}
}
close F;

sub betabinom {
	my ($sim, $a1, $b1, $a2, $b2) = @_;
	$a1 = $a1 + 1;
	$b1 = $b1 + 1;
	$a2 = $a2 + 1;
	$b2 = $b2 + 1;
	my @d1 = random_beta($sim, $a1, $b1);
	my @d2 = random_beta($sim, $a2, $b2);
	my $pos = 0;
	my $neg = 0;
	my $total = 0;
	my @dif = 0;
	for (my $i = 0; $i < @d1; $i++){
		my $subt = $d1[$i] - $d2[$i];
		push(@dif, $subt);
		if($subt > 0){
			$pos++;
		}
		elsif($subt < 0){
			$neg++;
		}
		$total++;
		if($i > 100){
			my $pp = ($pos)/$total;
			my $nn = ($neg)/$total;
			if($pp > 0.1 && $nn > 0.1){
				last;
			}
		}
	}
	my $obsp = ($pos)/$total;
	my $obsn = ($neg)/$total;
	my $pval;
	if($obsp < $obsn){
		$pval = (1 - $obsn + (1/$total));
	}
	elsif($obsn < $obsp){
		$pval = (1 - $obsp + (1/$total));
	}
	else{
		$pval = $obsp;
	}
	my $avedif = abs(mean(@dif));
	return($pval, $avedif);
}

sub mean {
	my (@vals) = @_;
	my $total = 0; 
	foreach(@vals){
		$total = $total + $_;
	}
	my $ave = $total/@vals;
	return($ave);
}
