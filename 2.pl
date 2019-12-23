#!/usr/bin/perl 
my $input = shift;
open FILE,$input;
open OUT,">network2.txt";
<FILE>;
while(<FILE>){
	chomp;
	my @a = split /\t/;
	$a[2]=~s/;.*//;
	$id{$a[0]}=$a[2];
	if($a[3]){
		if($a[1]=~/K\d+/){
			foreach my $t(split /; /,$a[3]){
				$ko{$t}{$a[0]}=1;
			}
			$o{$a[2]}="$a[2]\t$a[0]\t$a[1]\tProtein\n";
		}
		else{
			foreach my $t(split /; /,$a[3]){
				$co{$t}{$a[0]}=1;
			}
			$o{$a[2]}="$a[2]\t$a[0]\t$a[1]\tCompound\n";
		}
	}
}
print "KEGG pathway\tNo. of proteins\tNo. of compound\tProtein IDs\tCompound IDs\n";
foreach my $k(keys %ko){
	if($co{$k}){
		my @n1 = keys %{$ko{$k}};
		my @n2 = keys %{$co{$k}};
		foreach $a1(@n1){
			foreach $a2(@n2){
				$net{$a1}{$a2}++;
			}
		}
		my $n1=@n1;
		my $n2=@n2;
		print "$k\t$n1\t$n2\t@n1\t@n2\n";
	}
}
print OUT "node1\tnode2\tdegree\n";
foreach my $k1(keys %net){
	foreach my $k2(keys %{$net{$k1}}){
		print OUT "$id{$k1}\t$id{$k2}\t$net{$k1}{$k2}\n";
	}
}
close OUT;
open OUT,">table2.txt";
print OUT "Node name\tID\tKEGG ID\tType\n";
foreach my $k(keys %o){
	print OUT $o{$k};
}
