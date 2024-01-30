#!/perl -w
=head1 Name

	5.function_pre.pl -- A pipeline to annotation prokaryotic genome

=head1 Description

	This is a pipeline to annotation the prokaryote genome. There are several databases:
	KEGG,cazyDB,KEGG_2019,NOG,eggNOG,NR
	
=head1 Version

	Version: 1; by Jian-Yu Jiao

=head1 Usage

	perl 5.function_pre_V1.pl  [options] genome_path
		--KEGG			 
		--cazyDB
		--KEGG_2019
		--NOG
		--eggNOG
		--NR
		--cpu <int>         set the cpu number to use in parallel, default 20 for qsub and 5 for multi
		--help              output help information to screen
		
=head1 Example

	perl 5.function_pre_V1.pl  --KEGG_2019 xxxx_path
	perl 5.function_pre_V1.pl  -KEGG -KEGG_2019 -NOG -eggNOG -NR -cpu 20 /mdata/jiaojy/project/actino2/1.genomes_test_bins/remove_genome
=cut


use warnings;
use strict;

use Getopt::Long;



my($KEGG,$KEGG_2019,$NOG,$eggNOG,$NR,$cazyDB,$cpu,$Help);

GetOptions(
	"KEGG"=>\$KEGG,
	"cazyDB"=>\$cazyDB,
    "KEGG_2019"=>\$KEGG_2019,
	"NOG"=>\$NOG,
	"eggNOG"=>\$eggNOG,
	"NR"=>\$NR,
	"cpu:i"=>\$cpu,
	"help"=>\$Help
);
##get the introduction information
die `pod2text $0` if ($Help || @ARGV==0);

$cpu ||= 10;
my $path = shift;


opendir PATH, $path or die "can't open path $path";


foreach my $file (readdir PATH){
	if ($file =~ m/\.fna/ig){
	my @names = split/\.fna/, $file;
	my $folder = shift@names;
	system("mkdir $path/$folder/3.function_pre");	
                my $orf_new_pep = $folder . ".pep";
		my $orf_new_cds = $folder . ".cds";
		my $orf_cds_kegg =  $folder . ".kegg.daa";
		my $orf_cds_kegg_table =  $folder . ".kegg.table.out";
		my $orf_cds_kegg_table_parser =  $folder . ".kegg.table.parser.out";
		my $orf_cds_kegg_table_parser_rm_repeats =  $folder . ".kegg.table.parser.rm.repeats.out";
		my $orf_cds_kegg_table_parser_rm_repeats_KO =  $folder . ".kegg.table.parser.rm.repeats.KO.out";
		my $orf_cds_kegg_table_parser_rm_repeats_KO_level2 =  $folder . ".kegg.table.parser.rm.repeats.KO.level2.out";
		my $orf_cds_kegg_table_parser_rm_repeats_KO_level3 =  $folder . ".kegg.table.parser.rm.repeats.KO.level3.out";
		my $orf_cds_kegg_table_parser_rm_repeats_KO_level4 =  $folder . ".kegg.table.parser.rm.repeats.KO.level4.out";
		
		my $orf_cds_nog =  $folder . ".nog.daa";
		my $orf_cds_nog_table =  $folder . ".nog.table.out";
		my $orf_cds_nog_table_parser =  $folder . ".nog.table.parser.out";
		my $orf_cds_nog_table_parser_rm_repeats =  $folder . ".nog.table.parser.rm.repeats.out";
		my $orf_cds_nog_table_parser_rm_repeats_COG =  $folder . ".nog.table.parser.rm.repeats.COG.out";
		my $orf_cds_nog_table_parser_rm_repeats_COG_category =  $folder . ".nog.table.parser.rm.repeats.COG.category.out";
		my $orf_cds_nog_table_parser_rm_repeats_COG_category_category_sta =  $folder . ".nog.table.parser.rm.repeats.COG.category.sta";
		my $orf_cds_nog_table_parser_rm_repeats_COG_category_catalogue_sta =  $folder . ".nog.table.parser.rm.repeats.COG.catalogue.sta";
		
		my $orf_cds_nr =  $folder . ".nr.daa";
		my $orf_cds_nr_table =  $folder . ".nr.table.out";
		my $orf_cds_nr_table_parser =  $folder . ".nr.table.parser.out";
		my $orf_cds_nr_table_parser_rm_repeats =  $folder . ".nr.table.parser.rm.repeats.out";
		my $orf_cds_nr_table_parser_rm_repeats_taxon =  $folder . ".nr.table.parser.rm.repeats.taxon.out";
		
		my $cazy_out_pre =  $folder . ".";

#KEGG
if(defined $KEGG){
system("mkdir $path/$folder/3.function_pre/KEGG");
chdir("$path/$folder/3.function_pre/KEGG");	
	system("diamond_v0.7.9 blastx -d /software/database/kegg_136/kegg/genes.pep -q $path/$folder/2.gene_pre/$orf_new_cds -a $orf_cds_kegg -p $cpu -k 5 -e 1e-5");
    system("diamond_v0.7.9 view -a $orf_cds_kegg -o $orf_cds_kegg_table -f tab --compress 0");
    system("perl /mdata/jiaojy/script/kegg.parser/1.parse_diamond_table.pl $orf_cds_kegg_table /mdata/jiaojy/script/kegg.parser/1.genes.pep.title.parser.txt $orf_cds_kegg_table_parser");
    system("perl /mdata/jiaojy/script/kegg.parser/2.remove_repeats_in_m8_parser.pl $orf_cds_kegg_table_parser $orf_cds_kegg_table_parser_rm_repeats");
    system("perl /mdata/jiaojy/script/kegg.parser/3.get_ko_from_kegg_m8_table.pl $orf_cds_kegg_table_parser_rm_repeats /mdata/jiaojy/script/kegg.parser/3.ko_species2_ko_id.lian.txt $orf_cds_kegg_table_parser_rm_repeats_KO");
    system("perl /mdata/jiaojy/script/kegg.parser/4.get_level2_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO /mdata/jiaojy/script/kegg.parser/4.ko_parser.txt $orf_cds_kegg_table_parser_rm_repeats_KO_level2");
    system("perl /mdata/jiaojy/script/kegg.parser/5.get_level3_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO /mdata/jiaojy/script/kegg.parser/5.ko_parser.txt $orf_cds_kegg_table_parser_rm_repeats_KO_level3");
    system("perl /mdata/jiaojy/script/kegg.parser/6.get_level4_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO $orf_cds_kegg_table_parser_rm_repeats_KO_level4");
}

#KEGG-2019
#diamond v0.9
if(defined $KEGG_2019){
system("mkdir $path/$folder/3.function_pre/KEGG_2019");
chdir("$path/$folder/3.function_pre/KEGG_2019");	
	system("diamond_v0.9 blastx -d /software/database/KEGG_2019/kegg_2019.dmnd -q $path/$folder/2.gene_pre/$orf_new_cds -a $orf_cds_kegg -p $cpu -k 5 -e 1e-5");
    system("diamond_v0.9 view -a $orf_cds_kegg -o $orf_cds_kegg_table -f tab --compress 0");
    system("perl /mdata/jiaojy/script/kegg.parser/1.parse_diamond_table.pl $orf_cds_kegg_table /software/database/KEGG_2019/Metagenomes_Jun_2019.txt $orf_cds_kegg_table_parser");
    system("perl /mdata/jiaojy/script/kegg.parser/2.remove_repeats_in_m8_parser.pl $orf_cds_kegg_table_parser $orf_cds_kegg_table_parser_rm_repeats");
    system("perl /mdata/jiaojy/script/kegg.parser/3.get_ko_from_kegg_m8_table.pl $orf_cds_kegg_table_parser_rm_repeats /software/database/KEGG_2019/Metagenomes_Jun_2019_gene_ko_rm_em.txt $orf_cds_kegg_table_parser_rm_repeats_KO");
    system("perl /mdata/jiaojy/script/kegg.parser/4.get_level2_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO /mdata/jiaojy/script/kegg.parser/4.ko_parser.txt $orf_cds_kegg_table_parser_rm_repeats_KO_level2");
    system("perl /mdata/jiaojy/script/kegg.parser/5.get_level3_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO /mdata/jiaojy/script/kegg.parser/5.ko_parser.txt $orf_cds_kegg_table_parser_rm_repeats_KO_level3");
    system("perl /mdata/jiaojy/script/kegg.parser/6.get_level4_from_kegg_table.pl $orf_cds_kegg_table_parser_rm_repeats_KO $orf_cds_kegg_table_parser_rm_repeats_KO_level4");
}

#NOG
if(defined $NOG){
system("mkdir $path/$folder/3.function_pre/NOG");
chdir("$path/$folder/3.function_pre/NOG");	
    system("diamond_v0.7.9 blastx -d /software/database/nog_STRING_136/nog_STRING/protein.sequences.v9.05 -q $path/$folder/2.gene_pre/$orf_new_cds -a $orf_cds_nog -p $cpu -k 5 -e 1e-5");
	system("diamond_v0.7.9 view -a $orf_cds_nog -o $orf_cds_nog_table -f tab --compress 0");
    system("perl /mdata/jiaojy/script/nog.parser/1.parse_diamond_table.pl $orf_cds_nog_table /mdata/jiaojy/script/nog.parser/1.protein.sequences.v9.05.title.parser.txt $orf_cds_nog_table_parser");
    system("perl /mdata/jiaojy/script/nog.parser/2.remove_repeats_in_m8_parser.pl $orf_cds_nog_table_parser $orf_cds_nog_table_parser_rm_repeats");
    system("perl /mdata/jiaojy/script/nog.parser/3.get_cog_from_nog_m8_table.pl $orf_cds_nog_table_parser_rm_repeats /mdata/jiaojy/script/nog.parser/3.COG.mappings.v9.05.txt $orf_cds_nog_table_parser_rm_repeats_COG");
    system("perl /mdata/jiaojy/script/nog.parser/4.get_cog_category_info_from_cog_in_nog_table.pl $orf_cds_nog_table_parser_rm_repeats_COG /mdata/jiaojy/script/nog.parser/4.COG.funccat.txt $orf_cds_nog_table_parser_rm_repeats_COG_category");
    system("perl /mdata/jiaojy/script/nog.parser/5.COG_category_sta.pl $orf_cds_nog_table_parser_rm_repeats_COG_category $orf_cds_nog_table_parser_rm_repeats_COG_category_category_sta");
    system("perl /mdata/jiaojy/script/nog.parser/6.COG_catalogue_sta.pl $orf_cds_nog_table_parser_rm_repeats_COG_category $orf_cds_nog_table_parser_rm_repeats_COG_category_catalogue_sta");
}

#eggNOG
if(defined $eggNOG){
system("mkdir $path/$folder/3.function_pre/eggNOG");	
chdir("$path/$folder/3.function_pre/eggNOG");
    system("emapper.py -m diamond --cpu $cpu -i $path/$folder/2.gene_pre/$orf_new_pep -o $folder --seed_ortholog_evalue 0.00001 --temp_dir /mdata/jiaojy/tempDIR");
}

#NR
if(defined $NR){
system("mkdir $path/$folder/3.function_pre/NR");
chdir("$path/$folder/3.function_pre/NR");
	system("diamond_v0.7.9 blastx -d /software/database/nr_136/nr -q $path/$folder/2.gene_pre/$orf_new_cds -a $orf_cds_nr -p $cpu -k 5 -e 1e-5");
	system("diamond_v0.7.9 view -a $orf_cds_nr -o $orf_cds_nr_table -f tab --compress 0");
    system("perl /mdata/jiaojy/script/nr.parser/1.parse_diamond_table.pl $orf_cds_nr_table /mdata/jiaojy/script/nr.parser/1.nr.title.parser.txt $orf_cds_nr_table_parser");
    system("perl /mdata/jiaojy/script/nr.parser/2.remove_repeats_in_m8_parser.pl $orf_cds_nr_table_parser $orf_cds_nr_table_parser_rm_repeats");
    system("perl /mdata/jiaojy/script/nr.parser/3.get_taxonomy_for_blast_nr_by_gi.pl $orf_cds_nr_table_parser_rm_repeats /mdata/jiaojy/script/nr.parser/3.gi_taxid_prot.dmp /mdata/jiaojy/script/nr.parser/3.taxon_standard.txt $orf_cds_nr_table_parser_rm_repeats_taxon");
}

#cazyDB
if(defined $cazyDB){
system("mkdir $path/$folder/3.function_pre/CAZY");
chdir("$path/$folder/3.function_pre/");
	system("run_dbcan.py $path/$folder/2.gene_pre/$orf_new_pep protein --out_dir CAZY --out_pre $cazy_out_pre --db_dir /software/database/dbcan --dia_eval 1e-10 --dia_cpu $cpu --hmm_cpu $cpu --hotpep_cpu $cpu --use_signalP True");
}

	
}

}

