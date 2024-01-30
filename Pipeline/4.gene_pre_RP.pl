#!/usr/bin/perl

use warnings;
use strict;

opendir PATH, $ARGV[0] or die "can't open path $ARGV[0]";


my %hash = ();
my $path = $ARGV[0];

foreach my $file (readdir PATH){
	if ($file =~ m/\.fna/ig){
	my @names = split/\.fna/, $file;
	my $folder = shift@names;
	my $file_500 = join("_",$folder,"500.fna");
	system("mkdir $path/$folder/");
	system("mkdir $path/$folder/1.seq");
	system("mkdir $path/$folder/2.gene_pre");
	system("mkdir $path/$folder/2.gene_pre/1.marker_genes");
	system("cp $path/$file $path/$folder/1.seq");
	system("perl /software/Liworkshop/jiaojy_script/seqs.trim.pl $path/$folder/1.seq/$file $path/$folder/1.seq/$file_500 500");
	my $marker_stat = $folder . "_marker_gene.stat";
        my $orf_pep = $folder . ".raw.pep";
		my $orf_new_pep = $folder . ".pep";
        my $orf_cds = $folder . ".raw.cds";
		my $orf_new_cds = $folder . ".cds";
        my $orf_gff = $folder . ".raw.gff";
		my $orf_new_gff = $folder . ".gff";
        my $orf_stat = $folder . ".marker.raw.stat";
	system("cd $path/$folder/2.gene_pre");
	system("prodigal -a $path/$folder/2.gene_pre/$orf_pep -d $path/$folder/2.gene_pre/$orf_cds -f gff -g 11 -o $path/$folder/2.gene_pre/$orf_gff -p single -s /$path/$folder/2.gene_pre/$orf_stat -i $path/$folder/1.seq/$file_500");
		open (RAWPEP,"<$path/$folder/2.gene_pre/$orf_pep") or die "can't open $orf_pep";
		open (PEP, ">$path/$folder/2.gene_pre/$orf_new_pep") or die "can't open $orf_new_pep";
			my $seq_pep_num = 1;
			while(my $line = <RAWPEP>){
				$line =~ s/\s+$//ig;
				$line =~ s/\*//ig;
				if ($line =~ m/^>/ig){
					my $seq_pep_na = join("_",$folder,$seq_pep_num);
					print PEP ">$seq_pep_na\n";
					$seq_pep_num++;
				}else{print PEP "$line\n";}
			}
		close PEP;
		close RAWPEP;
		open (RAWCDS,"<$path/$folder/2.gene_pre/$orf_cds") or die "can't open $orf_cds";
                open (CDS, ">$path/$folder/2.gene_pre/$orf_new_cds") or die "can't open $orf_new_cds";
			my $seq_num = 1;
                        while(my $line = <RAWCDS>){
                                $line =~ s/\s+$//ig;
                                if ($line =~ m/^>/ig){
                                        my $seq_na = join("_",$folder,$seq_num);
                                        print CDS ">$seq_na\n";
                                        $seq_num++;
                                }else{print CDS "$line\n";}
                        }
		close CDS;
		close RAWCDS;
		
		#pre_marker_gene
			
		chdir("$path/$folder/2.gene_pre/");
		system("perl /software/cotton/AMPHORA2/Scripts/MarkerScanner_RP_jiao.pl  $path/$folder/2.gene_pre/$orf_new_pep");	
		system("mv $path/$folder/2.gene_pre/*.pep $path/$folder/2.gene_pre/1.marker_genes");
		system("mv $path/$folder/2.gene_pre/1.marker_genes/$orf_pep $path/$folder/2.gene_pre/");
		system("mv $path/$folder/2.gene_pre/1.marker_genes/$orf_new_pep $path/$folder/2.gene_pre/");
		system("cd $path/$folder/2.gene_pre/1.marker_genes");
		chdir("$path/$folder/2.gene_pre/1.marker_genes");
		system("seqkit stat -a *.pep > $orf_stat");
		chdir("$path");
		


}

}

