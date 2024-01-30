#!/usr/bin/perl
=head1 Name

	Phylogenetic_analysis_ISME_V1.1.pl -- A pipeline to Phylogenetic analysis

=head1 Description

	This is a pipeline to perform phylogenetic analysis
	      (1) brief genomic information
	      (2) check the quality of genomes
	      (3) reconstruct the Phylogenetic Tree based on 120 marker genes (GTDB-Tk and IQ-Tree)
	      (4) ANI analysis
	      (5) AAI analysis
	
=head1 Version

	Version: 1  ; by Jian-Yu Jiao; 2021.12.03
	Version: 1.1; by Jian-Yu Jiao; 2023.06.06

=head1 Usage

	perl Phylogenetic_analysis_ISME_V1.1.pl  [options] genome_path
		--stat             genomic information (seqkit; Shen et al., PloS one, 2016).
		--CheckM           check the completeness and contamination of each genome (checkM; Parks et al., Genome research, 2015).
		--Tree             reconstruct the Phylogenetic Tree based on 120 marker genes (GTDB-Tk and IQ-Tree, Jiao et al., ISME, 2021)
		--ANI              ANI analysis (Hua et al., NC, 2018; Jiao et al., ISME, 2021)
		--AAI              AAI analysis (pyANI; Pritchard et al., AM, 2016)
		--cpu <int>        set the cpu number to use in parallel, default 10
		--help             output help information to screen
		
=head1 Example

    Phylogenetic_analysis_ISME_V1.1.pl  -stat   -CheckM   -Tree    -ANI  -AAI  -cpu 10   xxxx_path
or  Phylogenetic_analysis_ISME_V1.1.pl  -stat   -CheckM   -Tree    -ANI  -AAI  -cpu 10   xxxx_path

=head1 References

    (1) Jiao et al., (2021). Microbial dark matter coming to light: challenges and opportunities. National Science Review, 8(3), nwaa280.
    (2) Jiao et al., (2021). Insight into the function and evolution of the Wood–Ljungdahl pathway in Actinobacteria. The ISME Journal, 1-14.
	(3) Hua et al., (2018). Genomic inference of the metabolism and evolution of the archaeal phylum Aigarchaeota. Nature communications, 9(1), 1-11.
	(4) Shen et al., (2016). SeqKit: a cross-platform and ultrafast toolkit for FASTA/Q file manipulation. PloS one, 11(10), e0163962.
	(5) Parks et al., (2015). CheckM: assessing the quality of microbial genomes recovered from isolates, single cells, and metagenomes. Genome research, 25(7), 1043-1055.
	(6) Pritchard et al. (2016). Genomics and taxonomy in diagnostics for food security: soft-rotting enterobacterial plant pathogens. Analytical Methods, 8(1), 12-24.

	
=cut


use warnings;
use strict;

use Getopt::Long;



my($stat,$CheckM,$Tree,$ANI,$AAI,$cpu,$Help);

GetOptions(
	"stat"=>\$stat,
	"CheckM"=>\$CheckM,
    "Tree"=>\$Tree,
	"ANI"=>\$ANI,
	"AAI"=>\$AAI,
	"cpu:i"=>\$cpu,
	"help"=>\$Help
);
##get the introduction information
die `pod2text $0` if ($Help || @ARGV==0);

$cpu ||= 10;
my $path = shift;
my @pa = split/\//, $path;
print ("@pa\n");
my $genomes_dir = pop @pa;
print ("@pa\n");

print "$genomes_dir\n";

my $genomes_dir_gtdb = $genomes_dir . "_gtdb";
print "$genomes_dir_gtdb\n";


my $genomes_dir_AAI = $genomes_dir . "_AAI";


my $results_path= join("/", @pa);
print "$results_path\n";

my $genomes_dir_results = $genomes_dir . "_results";
system("mkdir $results_path/$genomes_dir_results");


opendir PATH, $path or die "can't open path $path";

system("mkdir $results_path/All_fasta_500");


foreach my $file (readdir PATH){
	if ($file =~ m/\.fna/ig){
	my @names = split/\.fna/, $file;
	my $folder = shift@names;
	my $file_500 = join("_",$folder,"500.fna");
	system("mkdir $path/$folder/");
	system("mkdir $path/$folder/1.seq");
	system("mkdir $path/$folder/2.gene_pre");
	system("cp $path/$file $path/$folder/1.seq");
	system("perl /software/Liworkshop/jiaojy_script/seqs.trim.pl $path/$folder/1.seq/$file $path/$folder/1.seq/$file_500 500");
	system("cp $path/$folder/1.seq/$file_500 $results_path/All_fasta_500/$file");
	system("seqkit stat -a $path/$folder/1.seq/$file_500 > $path/$folder/1.seq/genome.stat"); 
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
	}
}

	

############################################################################################################################################
#stat
if(defined $stat){
chdir("$results_path/All_fasta_500");
system('seqkit stat -a *.fna > All.genomes.500.stat');
system("mkdir $results_path/$genomes_dir_results/01.stat");
system("mv All.genomes.500.stat $results_path/$genomes_dir_results/01.stat");
}
############################################################################################################################################

if(defined $CheckM){
chdir("$results_path");	
	system("checkm lineage_wf -f ./01.genomes_checkM.txt -t $cpu -x fna --tmpdir /mdata/Li_workshop/TEMP All_fasta_500 ./01.genomes_checkM");
    system("checkm qa -o 2 -f 01.genomes_checkM_GC.txt -t $cpu --tmpdir /mdata/Li_workshop/TEMP ./01.genomes_checkM/lineage.ms ./01.genomes_checkM");
    system("rm -r 01.genomes_checkM");
    system("rm -r 01.genomes_checkM.txt");
    system("mkdir $results_path/$genomes_dir_results/02.CheckM");
    system("mv 01.genomes_checkM_GC.txt $results_path/$genomes_dir_results/02.CheckM/checkM_results.txt");
}

############################################################################################################################################

if(defined $Tree){
chdir("$results_path");	
    system("gtdbtk classify_wf --skip_ani_screen --cpu $cpu --genome_dir All_fasta_500 --out_dir $genomes_dir_gtdb");
chdir("$results_path/$genomes_dir_gtdb/align");
	system("gunzip gtdbtk.bac120.user_msa.fasta.gz");
    system("iqtree2 -s gtdbtk.bac120.user_msa.fasta -alrt 1000 -B 1000 -T $cpu");
    system("mkdir $results_path/$genomes_dir_results/03.Tree");  
	system("mv $results_path/$genomes_dir_gtdb/align/gtdbtk.bac120.user_msa* $results_path/$genomes_dir_results/03.Tree");

}
############################################################################################################################################
if(defined $ANI){
chdir("$results_path");		
    system("average_nucleotide_identity.py -i All_fasta_500 -o All_fasta_500_ANI -g --gformat pdf,png -m ANIb --workers $cpu");
    system("mv All_fasta_500_ANI $results_path/$genomes_dir_results/04.ANI");

}
############################################################################################################################################
if(defined $AAI){
chdir("$results_path");		
    system("perl /software/Liworkshop/jiaojy_script/8.1.AAI_orthologs_calculate.pl $path");
    system("mkdir $results_path/$genomes_dir_results/05.AAI");     
    system("mv $genomes_dir_AAI/aai.txt $results_path/$genomes_dir_results/05.AAI");
    system("mv $genomes_dir_AAI/aai.matrix.table $results_path/$genomes_dir_results/05.AAI");
}

############################################################################################################################################
chdir("$results_path");		

#system("rm -r All_fasta_500");
#system("rm -r $genomes_dir_AAI");
#system("rm -r $genomes_dir_gtdb");
############################################################################################################################################
