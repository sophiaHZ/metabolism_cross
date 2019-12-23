#!/usr/bin/perl -w
if(@ARGV<5){
	print STDERR "perl $0 pathway quant ko enrichment tag [pathway_detail]\n";
	exit;
}
my $pathway = shift;
my $quant = shift;
my $ko = shift;
my $enrichment = shift;
my $tag = shift;
#my $pvalue = shift;
open OUT,">kegg_pathway_information.xls";

open FILE,$quant;
while(<FILE>){
	chomp;
	my @a = split /\t/;
	if($a[3]){
		$protein_quant{$a[0]}=$a[1];
		$protein_signl{$a[0]}=$a[3];
	}
}
open FILE,$ko;
while(<FILE>){
	chomp;
	my @a = split;
	if($a[1]){
		if($protein_signl{$a[0]}){
			$ko{$a[1]}{$a[0]}=1;
		}
	}
}
my %ko_signal;
foreach my $ko(keys %ko){
	my $prot_qu;
	my $prot_si;
	foreach my $pro(keys %{$ko{$ko}}){
		$prot_si.=$protein_signl{$pro};
		$ko_line{$ko}=1;
	}
	if($prot_si=~/up/i and $prot_si=~/down/i){
		$ko_signal{$ko}="yellow";
	}
	elsif($prot_si=~/down/i){
		$ko_signal{$ko}="green";
	}
	elsif($prot_si=~/up/i){
		$ko_signal{$ko}="red";
	}
}
open FILE,$enrichment;
<FILE>;
while(<FILE>){
	chomp;
	my @a = split /\t/;
	#if($a[0] eq 'Terms'){
	#	next;
	#}
	#if($a[6]<$pvalue){
	my @b = split /\s+/,$a[0];
	$enrich{$b[0]}=$_;
	#}
}

open FILE,$pathway;
while(<FILE>){
	chomp;
	if(/($tag\d+)/ and $enrich{$1}){
		my $pathway_id = $1;
		my $tt=$1;
		$tt=~s/$tag/map/;
		my $signal=0;
		while(my $my_line=<FILE>){
			chomp $my_line;
			unless($my_line){
				if($signal==0){
					$signal=1;
					next;
				}
				else{
					last;
				}
			}
			if($my_line=~/ko\:(K\d+)/ or $my_line=~/cpd\:(C\d+)/){
				if($ko_line{$1}){
					foreach my $id(keys %{$ko{$1}}){
						print OUT "$tt\t$id\t$1\t$ko_signal{$1}\t$protein_quant{$id}\t$protein_signl{$id}\n";
					}
				}
			}
			else{
				last;
			}

		}
	}
}
`python3 /home/swf/bin/metabolism/KEGG_color_map_new.py -i kegg_pathway_information.xls`;
