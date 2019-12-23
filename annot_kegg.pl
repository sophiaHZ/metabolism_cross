#!/usr/bin/perl -w
use strict;
if(@ARGV<2){
	print STDERR "perl $0 pathway [ko_list] tag\n";
	exit;
}
my $pathway = shift;
my $ko_list = shift;
my $tag = shift;
open FILE,$pathway;
my $name;
my %name_id;
my %hash;
my %gene_desc;
my %level1;
my %level2;
while(<FILE>){
	chomp;
	if(/^($tag\d+ (.*)) \(\d+\)/){
		$name = $1;
#		$name=~s/ - .*//;
	}
	elsif(/ko:(K\d+) (.*)/){
		$hash{$1}{$name}=1;
		$gene_desc{$1}=$2;
	}
	elsif(/cpd:(C\d+) (.*)/){
		$hash{$1}{$name}=1;
		$gene_desc{$1}=$2;
	}
}
open FILE,$ko_list;
print "Protein accession\tKEGG KO No.\tKEGG Gene\tKEGG pathway\n";
while(<FILE>){
	chomp;
	print;
	my @a = split /\t/;
	if($a[1] and $hash{$a[1]}){
		print "\t$gene_desc{$a[1]}\t";
		foreach my $k(keys %{$hash{$a[1]}}){
			print "$k; ";
		}
	}
	print "\n";
}
