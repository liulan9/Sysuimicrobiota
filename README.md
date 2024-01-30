# Sysuimicrobiota
# Function annotation
## use
nohup perl 5.function_pre_V1.pl -cazyDB -KEGG_2019 -NOG -eggNOG -NR -cpu 20 genomes &

# Tree-GTDB
## use
nohup Phylogenetic_analysis_ISME_V1.1.pl  -Tree -cpu 15 genomes &

# Tree-Sixteen ribosomal protein sequences
## use
nohup perl gene_pre_RP.pl genomes &
nohup perl get_all_marker_gene_RP.pl genomes &
cd genomes_marker_genes
cd all_marker_gene/
catfasta2phyml.pl -f --sequential --concatenate *.trimal.fa > genomes.all.16rps.fasta
nohup iqtree_1.6.10 -s genomes.all.16rps.fasta -alrt 1000 -bb 1000 -nt 10 &

# Tree-31 conserved marker genes
## use
nohup perl gene_pre.pl genomes & 
nohup perl get_all_marker_gene.pl genomes &
cd genomes_marker_genes
cd all_marker_gene/
trimal -in dnaG.pep.align -out dnaG.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout dnaG.pep.align.trimal.html
trimal -in frr.pep.align -out frr.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout frr.pep.align.trimal.html
trimal -in infC.pep.align -out infC.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout infC.pep.align.trimal.html
trimal -in nusA.pep.align -out nusA.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout nusA.pep.align.trimal.html
trimal -in pgk.pep.align -out pgk.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout pgk.pep.align.trimal.html
trimal -in pyrG.pep.align -out pyrG.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout pyrG.pep.align.trimal.html
trimal -in rplA.pep.align -out rplA.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplA.pep.align.trimal.html
trimal -in rplB.pep.align -out rplB.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplB.pep.align.trimal.html
trimal -in rplC.pep.align -out rplC.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplC.pep.align.trimal.html
trimal -in rplD.pep.align -out rplD.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplD.pep.align.trimal.html
trimal -in rplE.pep.align -out rplE.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplE.pep.align.trimal.html
trimal -in rplF.pep.align -out rplF.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplF.pep.align.trimal.html
trimal -in rplK.pep.align -out rplK.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplK.pep.align.trimal.html
trimal -in rplL.pep.align -out rplL.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplL.pep.align.trimal.html
trimal -in rplM.pep.align -out rplM.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplM.pep.align.trimal.html
trimal -in rplN.pep.align -out rplN.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplN.pep.align.trimal.html
trimal -in rplP.pep.align -out rplP.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplP.pep.align.trimal.html
trimal -in rplS.pep.align -out rplS.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplS.pep.align.trimal.html
trimal -in rplT.pep.align -out rplT.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rplT.pep.align.trimal.html
trimal -in rpmA.pep.align -out rpmA.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpmA.pep.align.trimal.html
trimal -in rpoB.pep.align -out rpoB.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpoB.pep.align.trimal.html
trimal -in rpsB.pep.align -out rpsB.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsB.pep.align.trimal.html
trimal -in rpsC.pep.align -out rpsC.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsC.pep.align.trimal.html
trimal -in rpsE.pep.align -out rpsE.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsE.pep.align.trimal.html
trimal -in rpsI.pep.align -out rpsI.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsI.pep.align.trimal.html
trimal -in rpsJ.pep.align -out rpsJ.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsJ.pep.align.trimal.html
trimal -in rpsK.pep.align -out rpsK.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsK.pep.align.trimal.html
trimal -in rpsM.pep.align -out rpsM.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsM.pep.align.trimal.html
trimal -in rpsS.pep.align -out rpsS.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout rpsS.pep.align.trimal.html
trimal -in smpB.pep.align -out smpB.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout smpB.pep.align.trimal.html
trimal -in tsf.pep.align -out tsf.pep.align.trimal.fa -gt 0.95 -cons 50 -htmlout tsf.pep.align.trimal.html
catfasta2phyml.pl -f -sequential -concatenate *.trimal.fa > genomes.all.maker.fasta
nohup iqtree_1.6.10 -s genomes.all.maker.fasta -alrt 1000 -bb 1000 -nt 10 &
