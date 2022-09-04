#About: This script analyses RNAseq reads
#Author: Lilian Mbaisi
#Date created: June 15th, 2022
#Date last modified: July 8th,2022
#Prerequisites:
	#pycoQC
	#porechop
	#sortmerna
	#isonclust
	#isoncorrect
	#minimap
	#nanocount
#RNA‐seq data analyses typically consist of :
	#(1) accurate mapping of millions of short sequencing reads to a reference genome, including the identification of splicing events; 
	#(2) quantifying expression levels of genes, transcripts, and exons; 
	#(3) differential analysis of gene expression among different biological conditions; and 
	#(4) biological interpretation of differentially expressed genes. 
#######################################################################################
#STEP ONE: PYCOQC
###################

#######################################################################################
#STEP TWO: PORECHOP
###################
#Citing porechop:Wick, R. Porechop: adapter trimmer for Oxford Nanopore reads 2018 https://github.com/rrwick/Porechop
################################
#Combo script that runs porechop on all fastq files and concatenates them into one porechopped fastq file per barcode
for n in {01..12};do mkdir output_reads2/barcode$n;(for i in FB_RNAseq_passed/barcode$n/*.fastq.gz;do base=$(basename $i -s="FB_RNAseq_passed/barcode*/*.fastq.gz");porechop -i $i --format fastq.gz -o output_reads2/barcode$n/$base.porechopped.fastq; cat output_reads2/barcode$n/$base.porechopped.fastq >> output_reads2/barcode$n/barcode$n.porechopped.fastq;done);done

#concatenate output files per barcode
for n in {01..12};do mkdir output_reads/barcode$n;(for i in FB_RNAseq_passed/barcode$n/*.fastq.gz;do base=$(basename $i -s="FB_RNAseq_passed/barcode*/*.fastq.gz"); cat output_reads2/barcode$n/$base.porechopped.fastq >> output_reads/barcode$n/barcode$n.porechopped.fastq;done);done

#######################################################################################
#STEP3: SORTMERNA
###################
#Citing sortmerna:Kopylova E., Noé L. and Touzet H., "SortMeRNA: Fast and accurate filtering of ribosomal RNAs in metatranscriptomic data", Bioinformatics (2012), doi: 10.1093/bioinformatics/bts611. 
################################
#export PATH=/home/lily/Desktop/Jackie_New_RNASeq_Data_June15/bin:$PATH
#Install sortmerna version >3

for n in {01..12};do sortmerna --ref sortmerna/data/rRNA_databases/silva-arc-16s-id95.fasta --ref sortmerna/data/rRNA_databases/silva-arc-23s-id98.fasta --ref sortmerna/data/rRNA_databases/silva-bac-16s-id90.fasta --ref sortmerna/data/rRNA_databases/silva-bac-23s-id98.fasta --ref sortmerna/data/rRNA_databases/silva-euk-18s-id95.fasta --ref sortmerna/data/rRNA_databases/silva-euk-28s-id98.fasta -reads final_outputs/barcode$i.porechopped.fastq --aligned aligned_barcode$i.porechopped --fastx --other unaligned_barcode$i.porechopped --blast '1 cigar qcov' barcode$i.blast.out;rm -rf /home/lily/sortmerna/run/;done

#######################################################################################
#STEP 4: ISONCLUST & ISONCORRECT
################################
#isonclust citation: Kristoffer Sahlin, Paul Medvedev. De Novo Clustering of Long-Read Transcriptome Data Using a Greedy, Quality-Value Based Algorithm, Journal of Computational Biology 2020, 27:4, 472-484

#isoncorrect citation: Sahlin, K., Medvedev, P. Error correction enables use of Oxford Nanopore technology for reference-free transcriptome analysis. Nat Commun 12, 2 (2021). https://doi.org/10.1038/s41467-020-20340-8 Link.
################################

#Should be run on unzipped fastq
for i in *.porechopped.fq; do base=$(basename $i .porechopped.fq); ./isONclust write_fastq --clusters $base.isonclust.output/final_clusters.tsv --fastq $i --outfolder $base.isonclust.fastq.output --N 1;isONcorrect/./run_isoncorrect --fastq_folder $base.isonclust.fastq.output --k 13 --w 20 --outfolder $base.correction;done
#move isoncorrect files to isonclust before running above.

# OPTIONAL BELOW TO MERGE ALL CORRECTED READS INTO ONE FILE
for outfolder in *.correction; do base=$(basename $outfolder .correction);
touch $outfolder/all_corrected_reads.fq

OUTFILES=$outfolder/*/corrected_reads.fastq

for f in $OUTFILES
do 
  echo $f
  cat $f >> $outfolder/$base.all_corrected_reads.fq
done

echo
echo "Finished with pipeline and wrote corrected reads to: " $outfolder/$base.all_corrected_reads.fq
echo ;

done

#copy corrected reads to one folder
mkdir all_corrected_reads

for outfolder in *.correction; do base=$(basename $outfolder .correction); cp $outfolder/$base.all_corrected_reads.fq all_corrected_reads/;echo "Finished copying corrected reads" $base.all_corrected_reads.fq "to: all_corrected_reads/";done
#######################################################################################
#STEP 5 MINIMAP2 and NANOCOUNT: Mapping to reference and countinh
######################################
#Citing minimap2:
    #Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34:3094-3100. doi:10.1093/bioinformatics/bty191

    #Li, H. (2021). New strategies to improve minimap2 alignment accuracy. Bioinformatics, 37:4572-4574. doi:10.1093/bioinformatics/btab705
#Citing nanocount:
	#Josie Gleeson, Adrien Leger, Yair D J Prawer, Tracy A Lane, Paul J Harrison, Wilfried Haerty, Michael B Clark, Accurate expression quantification from nanopore direct RNA sequencing with NanoCount, Nucleic Acids Research, 2021;, gkab1129, https://doi.org/10.1093/nar/gkab1129
################################
#Mapping to host genome and counting
for i in *.all_corrected_reads.fq; do minimap2 -a -x splice -k13 AR.fa $i > output_reads/$i.aln.sam;done
for i in *.aln.sam; do base=$(basename $i .aln.sam);samtools view -bS $i > $base.bam;done
for i in *.bam; do base=$(basename $i .bam); samtools view -h -F 4 $i -o $base.mapped.bam; samtools sort $base.mapped.bam > $base.mapped.sorted.bam; samtools index $base.mapped.sorted.bam;done
for i in *.mapped.sorted.bam; do base=$(basename $i .mapped.sorted.bam); NanoCount -i $i -o $base.nanocount.tsv -b $base.nanocount.bam -t 0.8 --verbose;done

#Mapping to host transcriptome and counting

for i in *.all_corrected_reads.fq; do base=$(basename $i .all_corrected_reads.fq);minimap2 -a -x map-ont -k13 GCF_016920715.1_AaraD3_rna.fna $i > output_reads_rna/$base.aln.rna.sam;done
for i in *.aln.rna.sam; do base=$(basename $i .aln.rna.sam);samtools view -bS $i > $base.bam;done
for i in *.bam; do base=$(basename $i .bam); samtools view -h -F 4 $i -o $base.mapped.bam; samtools sort $base.mapped.bam > $base.mapped.sorted.bam; samtools index $base.mapped.sorted.bam;done
for i in *.mapped.sorted.bam; do base=$(basename $i .mapped.sorted.bam); NanoCount -i $i -o $base.nanocount.tsv -b $base.nanocount.bam -t 0.8 --verbose --extra_tx_info;done


#######################################################################################
#HTSEQ-COUNT JULY 1ST- Work in progress
for i in *.bam; do echo $i; htseq-count $i --idattr=gbkey ~/Desktop/RNA_Seq_work/GCF_016920715.1_AaraD3_genomic.gff;done
########################################################################################





