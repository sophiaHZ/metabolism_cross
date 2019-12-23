#!/usr/bin/perl -w
use strict;
use Text::NSP::Measures::2D::Fisher::right;
use Statistics::ChisqIndep;
use Excel::Writer::XLSX;
use Statistics::R;
if(@ARGV<2){
	print STDERR "perl $0 pathway_annotation dif_protein_list\n";
	exit;
}
my $pathway_annotation = shift;
my $dif_protein_list = shift;
my $signi = shift;
my $min_number = shift;
my $type = shift;
my $pathway_level_file="/home/pxj/database/KEGG/all_pathway_level.txt";
my %pathway_level;
open FILE,$pathway_level_file;
while(<FILE>){
	chomp;
	my @a = split /\t/;
	if($a[3] eq '1.0 Global and overview maps'){
		$a[0]=~s/map//;
		$pathway_level{$a[0]}=$a[3];
	}
}
unless($signi){
	$signi=0;
}
$signi-=1;
unless($min_number){
	$min_number=1.5;
}
open FILE,$dif_protein_list;
my %sing;
my %quant;
my $dif_header=<FILE>;
chomp $dif_header;
my @dif_header=split /\t/,$dif_header;
my %quant1;
while(<FILE>){
	chomp;
	my @a  = split /\t/;
	$quant{$a[0]}=$_;
	if($a[3]){
		$quant1{$a[0]}{"all"}{$a[1]}=$_;
		$sing{$a[0]}{all}=1;
	}
}
open FILE,$pathway_annotation;
<FILE>;
my %name;
my %dif;
my %dif_n;
my %all_n;
my %all;
my $pathway_tag;
while(<FILE>){
	chomp;
	s/'s//g;
	my @a = split /\t/;
	if($quant{$a[0]} and $a[3]){
		my @b = split /; /,$a[3];
		my $k = $a[0];
		foreach my $pathway_tag(@b){
			$pathway_tag=~s/ \(\d+\)$//;
			$pathway_tag=~s/ - [^\-]*$//;
			$all_n{$pathway_tag}{$k}=1;
			$all{$k}=1;
			if($sing{$k}){
				foreach my $kn(keys %{$sing{$k}}){
					$dif_n{$kn}{$pathway_tag}{$k}=1;
					$dif{$kn}{$k}=1;
				}
			}
		}
	}
}
foreach my $k(keys %dif){
	KEGG($k,$type);
}
sub KEGG{
	my $sigl_tag = shift;
	my $file_tag = shift;
	my $chi = new Statistics::ChisqIndep;
	my %out_line;
	my %out;
	my %rank;
	my %fold;
	my $mm=0;
	my $mm1=0;
	my $mm2=0;
	my $npp=keys %all;
	my $np1=keys %{$dif{$sigl_tag}};
	open OUT,">${file_tag}_pathway_enrichment.xls";
#	print OUT "Type\tKEGG pathway\tMapping\tBackground\tAll Mapping\tAll Background\tFold enrichment\tp value\tRelated proteins\n";
	foreach my $k(keys %{$dif_n{$sigl_tag}}){
		my $n11 = keys %{$dif_n{$sigl_tag}{$k}};
		my @n11 = keys %{$dif_n{$sigl_tag}{$k}};
		my $n1p = keys %{$all_n{$k}};
		my $kd=$k;
		$kd=~s/[a-z]*//;
		$kd=~s/ .*//;
		if($pathway_level{$kd}){
			next;
		}
		else{
			my $p_value = calculateStatistic(n11=>$n11+$signi,n1p=>$n1p,np1=>$np1+$signi,npp=>$npp);
			if($p_value==0){
				my @load_data = ([$n11,$np1],[$n1p,$npp]);
				$chi->load_data(\@load_data);
				$p_value=$chi->{p_value};
			}
			my $fold = sprintf("%.2f",$n11/$n1p/$np1*$npp);
			my $log_fold = sprintf("%.2f",-log($p_value)/log(10));
			$out_line{$k}="$file_tag\t$k\t$n11\t$n1p\t$np1\t$npp\t$fold\t$p_value\t@n11";
#			$out{$k}="$k\t$n11\t$n1p\t$np1\t$npp\t$fold\t$p_value\t$log_fold\t@n11";
			$rank{$k}=$p_value;
#			$fold{$k}=$fold;
		}
	}
	foreach my $k(sort{$rank{$a}<=>$rank{$b}} keys %rank){
		print OUT "$out_line{$k}\n";
	}
	close OUT;
}
