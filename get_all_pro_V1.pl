#!/perl -w
=head1 Name

	4.3_get_all_pro_V1.pl -- get pep from results

=head1 Description

	Get pep sequences from results of KEGG,cazyDB,KEGG_2019,NOG,eggNOG,NR
	
=head1 Version

	Version: 1; by Jian-Yu Jiao

=head1 Usage

	perl 4.3_get_all_pro_V1.pl  [options] genome_path
		--KEGG			 
		--cazyDB
		--KEGG_2019
		--NOG
		--eggNOG
		--NR
		--help              output help information to screen
		
=head1 Example

	perl 4.3_get_all_pro_V1.pl  --KEGG_2019 K00001,K00002  xxxx_path
	perl 5.function_pre_V1.pl  -KEGG -KEGG_2019 -NOG -eggNOG -NR -cazyDB -cpu 20 /mdata/jiaojy/project/actino2/1.genomes_test_bins/remove_genome
=cut
use warnings;
use strict;
use Getopt::Long;

my($KEGG,$KEGG_2019,$NOG,$eggNOG,$NR,$cazyDB,$cpu,$Help);

GetOptions(
	"KEGG:s"=>\$KEGG,
	"cazyDB:s"=>\$cazyDB,
    "KEGG_2019:s"=>\$KEGG_2019,
	"NOG:s"=>\$NOG,
	"eggNOG:s"=>\$eggNOG,
	"NR:s"=>\$NR,
	"help"=>\$Help
);
##get the introduction information
die `pod2text $0` if ($Help || @ARGV==0);
my $path = shift;


opendir PATH, $path or die "can't open path $path";

my %hash = ();
my $new_path = join("_",$path,"pro");
system("mkdir $new_path");
system("mkdir $new_path/all_pep");
foreach my $file (readdir PATH){
	if ($file =~ m/\.fna/ig){
	my @names = split/\.fna/, $file;
	my $folder = shift@names;
	my $pep_path = join ("\/",$path,$folder,"2.gene_pre");
	my $pep_file = join ("\.",$folder,"pep");
		opendir PEPPATH, $pep_path or die "can't open path $path/$folder/2.gene_pre";
		foreach my $pep (readdir PEPPATH){
			if($pep =~ m/$pep_file/ig){
				system("cp $pep_path/$pep $new_path/all_pep");
			}
		}

		
if(defined $KEGG_2019){
	my $anno_path = join ("\/",$path,$folder,"3.function_pre","KEGG_2019");
	my $anno_file = join ("\.",$folder,"kegg.table.parser.rm.repeats.KO.out");
	my @kos = split/,/,$KEGG_2019;		
	open (KOFILE, "<$path/$folder/3.function_pre/KEGG_2019/$anno_file") or die "can't open file $anno_file";
		while (my $line = <KOFILE>){
			$line =~ s/\s+$//ig;
			my @list = split/\s+/,$line;
			foreach my $ko(@kos){
				if($list[6] =~ $ko){
					open (KOOUT, ">>$new_path/$ko") or die "can't open file $ko";
					print KOOUT "$list[0]\n";
				}
		
			}
		}	
}
		
	}
}

system("cat $new_path/all_pep/*.pep > $new_path/all.pep");

if(defined $KEGG_2019){
	chdir("$new_path");
	my @kos = split/,/,$KEGG_2019;	
	foreach my $ko(@kos){
		my $ko_fasta =  $ko . ".pep.fa";
		system("seqkit grep -f $ko all.pep > $ko_fasta");
	}
}
	


