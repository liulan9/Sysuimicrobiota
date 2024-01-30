#!/usr/bin/perl

use warnings;
use strict;

opendir PATH, $ARGV[0] or die "can't open path $ARGV[0]";

my %hash = ();
my $path = $ARGV[0];
my $new_path = join("_",$path,"marker_genes");
system("mkdir $new_path");


foreach my $file (readdir PATH){
	if ($file =~ m/\.fna/ig){
	my @names = split/\.fna/, $file;
	my $folder = shift@names;
	my $pep_path = join ("\/",$path,$folder,"2.gene_pre","1.marker_genes");
		opendir PEPPATH, $pep_path or die "can't open path $path/$folder/2.gene_pre/1.marker_genes";
		foreach my $pep (readdir PEPPATH){
			if($pep =~ m/pep/ig){
				my $new_pep = join("_",$folder,$pep);
					open (MARKER,"<$path/$folder/2.gene_pre/1.marker_genes/$pep");
					open (OUT,">$new_path/$new_pep");
						while(my $line = <MARKER>){
							$line =~ s/\s+$//ig;
							if ($line =~ s/^>//ig){
								my @list_line = split/_/,$line;
								my $remove = pop@list_line;
								my $new_name = join("_",@list_line);
								print OUT ">$new_name\n";
							}else {print OUT "$line\n";}
						}
			}
		}
	}	

}
system("mkdir $new_path/all_marker_gene");
chdir("$new_path/all_marker_gene");
system("cat ../*rpL14.pep > rpL14.pep");
system("cat ../*rpL15.pep > rpL15.pep");
system("cat ../*rpL16.pep > rpL16.pep");
system("cat ../*rpL18.pep > rpL18.pep");
system("cat ../*rpL2.pep > rpL2.pep");
system("cat ../*rpL22.pep > rpL22.pep");
system("cat ../*rpL24.pep > rpL24.pep");
system("cat ../*rpL3.pep > rpL3.pep");
system("cat ../*rpL4.pep > rpL4.pep");
system("cat ../*rpL5.pep > rpL5.pep");
system("cat ../*rpL6.pep > rpL6.pep");
system("cat ../*rpS10.pep > rpS10.pep");
system("cat ../*rpS17.pep > rpS17.pep");
system("cat ../*rpS19.pep > rpS19.pep");
system("cat ../*rpS3.pep > rpS3.pep");
system("cat ../*rpS8.pep > rpS8.pep");

system("muscle -maxiters 100 -in rpL14.pep -out rpL14.pep.align -quiet");
system("muscle -maxiters 100 -in rpL15.pep -out rpL15.pep.align -quiet");
system("muscle -maxiters 100 -in rpL16.pep -out rpL16.pep.align -quiet");
system("muscle -maxiters 100 -in rpL18.pep -out rpL18.pep.align -quiet");
system("muscle -maxiters 100 -in rpL2.pep -out rpL2.pep.align -quiet");
system("muscle -maxiters 100 -in rpL22.pep -out rpL22.pep.align -quiet");
system("muscle -maxiters 100 -in rpL24.pep -out rpL24.pep.align -quiet");
system("muscle -maxiters 100 -in rpL3.pep -out rpL3.pep.align -quiet");
system("muscle -maxiters 100 -in rpL4.pep -out rpL4.pep.align -quiet");
system("muscle -maxiters 100 -in rpL5.pep -out rpL5.pep.align -quiet");
system("muscle -maxiters 100 -in rpL6.pep -out rpL6.pep.align -quiet");
system("muscle -maxiters 100 -in rpS10.pep -out rpS10.pep.align -quiet");
system("muscle -maxiters 100 -in rpS17.pep -out rpS17.pep.align -quiet");
system("muscle -maxiters 100 -in rpS19.pep -out rpS19.pep.align -quiet");
system("muscle -maxiters 100 -in rpS3.pep -out rpS3.pep.align -quiet");
system("muscle -maxiters 100 -in rpS8.pep -out rpS8.pep.align -quiet");


system("trimal -in rpL14.pep.align -out rpL14.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL14.pep.align.trimal.html");
system("trimal -in rpL15.pep.align -out rpL15.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL15.pep.align.trimal.html");
system("trimal -in rpL16.pep.align -out rpL16.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL16.pep.align.trimal.html");
system("trimal -in rpL18.pep.align -out rpL18.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL18.pep.align.trimal.html");
system("trimal -in rpL2.pep.align -out rpL2.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL2.pep.align.trimal.html");
system("trimal -in rpL22.pep.align -out rpL22.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL22.pep.align.trimal.html");
system("trimal -in rpL24.pep.align -out rpL24.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL24.pep.align.trimal.html");
system("trimal -in rpL3.pep.align -out rpL3.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL3.pep.align.trimal.html");
system("trimal -in rpL4.pep.align -out rpL4.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL4.pep.align.trimal.html");
system("trimal -in rpL5.pep.align -out rpL5.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL5.pep.align.trimal.html");
system("trimal -in rpL6.pep.align -out rpL6.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpL6.pep.align.trimal.html");
system("trimal -in rpS10.pep.align -out rpS10.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpS10.pep.align.trimal.html");
system("trimal -in rpS17.pep.align -out rpS17.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpS17.pep.align.trimal.html");
system("trimal -in rpS19.pep.align -out rpS19.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpS19.pep.align.trimal.html");
system("trimal -in rpS3.pep.align -out rpS3.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpS3.pep.align.trimal.html");
system("trimal -in rpS8.pep.align -out rpS8.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpS8.pep.align.trimal.html");









