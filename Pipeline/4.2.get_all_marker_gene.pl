#!/usr/bin/perl
#
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
system("cat ../*dnaG.pep > dnaG.pep");
system("cat ../*frr.pep > frr.pep");
system("cat ../*infC.pep > infC.pep");
system("cat ../*nusA.pep > nusA.pep");
system("cat ../*pgk.pep > pgk.pep");
system("cat ../*pyrG.pep > pyrG.pep");
system("cat ../*rplA.pep > rplA.pep");
system("cat ../*rplB.pep > rplB.pep");
system("cat ../*rplC.pep > rplC.pep");
system("cat ../*rplD.pep > rplD.pep");
system("cat ../*rplE.pep > rplE.pep");
system("cat ../*rplF.pep > rplF.pep");
system("cat ../*rplK.pep > rplK.pep");
system("cat ../*rplL.pep > rplL.pep");
system("cat ../*rplM.pep > rplM.pep");
system("cat ../*rplN.pep > rplN.pep");
system("cat ../*rplP.pep > rplP.pep");
system("cat ../*rplS.pep > rplS.pep");
system("cat ../*rplT.pep > rplT.pep");
system("cat ../*rpmA.pep > rpmA.pep");
system("cat ../*rpoB.pep > rpoB.pep");
system("cat ../*rpsB.pep > rpsB.pep");
system("cat ../*rpsC.pep > rpsC.pep");
system("cat ../*rpsE.pep > rpsE.pep");
system("cat ../*rpsI.pep > rpsI.pep");
system("cat ../*rpsJ.pep > rpsJ.pep");
system("cat ../*rpsK.pep > rpsK.pep");
system("cat ../*rpsM.pep > rpsM.pep");
system("cat ../*rpsS.pep > rpsS.pep");
system("cat ../*smpB.pep > smpB.pep");
system("cat ../*tsf.pep > tsf.pep");

system("muscle -maxiters 100 -in dnaG.pep -out dnaG.pep.align");
system("muscle -maxiters 100 -in frr.pep -out frr.pep.align");
system("muscle -maxiters 100 -in infC.pep -out infC.pep.align");
system("muscle -maxiters 100 -in nusA.pep -out nusA.pep.align");
system("muscle -maxiters 100 -in pgk.pep -out pgk.pep.align");
system("muscle -maxiters 100 -in pyrG.pep -out pyrG.pep.align");
system("muscle -maxiters 100 -in rplA.pep -out rplA.pep.align");
system("muscle -maxiters 100 -in rplB.pep -out rplB.pep.align");
system("muscle -maxiters 100 -in rplC.pep -out rplC.pep.align");
system("muscle -maxiters 100 -in rplD.pep -out rplD.pep.align");
system("muscle -maxiters 100 -in rplE.pep -out rplE.pep.align");
system("muscle -maxiters 100 -in rplF.pep -out rplF.pep.align");
system("muscle -maxiters 100 -in rplK.pep -out rplK.pep.align");
system("muscle -maxiters 100 -in rplL.pep -out rplL.pep.align");
system("muscle -maxiters 100 -in rplM.pep -out rplM.pep.align");
system("muscle -maxiters 100 -in rplN.pep -out rplN.pep.align");
system("muscle -maxiters 100 -in rplP.pep -out rplP.pep.align");
system("muscle -maxiters 100 -in rplS.pep -out rplS.pep.align");
system("muscle -maxiters 100 -in rplT.pep -out rplT.pep.align");
system("muscle -maxiters 100 -in rpmA.pep -out rpmA.pep.align");
system("muscle -maxiters 100 -in rpoB.pep -out rpoB.pep.align");
system("muscle -maxiters 100 -in rpsB.pep -out rpsB.pep.align");
system("muscle -maxiters 100 -in rpsC.pep -out rpsC.pep.align");
system("muscle -maxiters 100 -in rpsE.pep -out rpsE.pep.align");
system("muscle -maxiters 100 -in rpsI.pep -out rpsI.pep.align");
system("muscle -maxiters 100 -in rpsJ.pep -out rpsJ.pep.align");
system("muscle -maxiters 100 -in rpsK.pep -out rpsK.pep.align");
system("muscle -maxiters 100 -in rpsM.pep -out rpsM.pep.align");
system("muscle -maxiters 100 -in rpsS.pep -out rpsS.pep.align");
system("muscle -maxiters 100 -in smpB.pep -out smpB.pep.align");
system("muscle -maxiters 100 -in tsf.pep -out tsf.pep.align");
system('seqkit stat -a *.pep.align > all.align.stat');

