#!/usr/bin/env bash

# #quality check

# FILE="/home/kebaso/Desktop/Mini_projects/nextflow/data/CK/barcode01/*.gz"

# #echo $(basename -a $FILE)
# for file in ${FILE}
# do
#         f="$(basename -a $file)"
#         `mkdir -p results`
#         nanoQC -o results ${file}
#         `mv results/nanoQC.html results/${f}.html`
        
        
# done

#Combo script that runs porechop on all fastq files and concatenates them into one porechopped fastq 


# for n in {02..03}
# do 
# mkdir -p output_reads2/barcode$n


# for i in CK/barcode$n/*.fastq.gz;
# do 
# #echo ${i}
# base=$(basename $i -s="FB_RNAseq_passed/barcode*/*.fastq.gz");
# #echo ${base}

# porechop -i $i --format fastq -o output_reads2/barcode$n/$base.porechopped.fastq


# cat output_reads2/barcode$n/$base.porechopped.fastq >> output_reads2/barcode$n/barcode$n.porechopped.fastq


# done
# done

# #move output files per barcode
# for n in {02..03};
# do 
# mkdir -p output_reads;
# mv output_reads2/barcode$n/barcode$n* ./output_reads/
# done

# # #STEP3: SORTMERNA

# # #SORTMERNA

# for i in {02..03}; do 
#     sortmerna \
#         --ref ./rRNA_databases/silva-arc-16s-id95.fasta \
#         --ref ./rRNA_databases/silva-arc-23s-id98.fasta \
#         --ref ./rRNA_databases/silva-bac-16s-id90.fasta \
#         --ref ./rRNA_databases/silva-bac-23s-id98.fasta \
#         --ref ./rRNA_databases/silva-euk-18s-id95.fasta \
#         --ref ./rRNA_databases/silva-euk-28s-id98.fasta \
#         -reads output_reads/barcode$i.porechopped.fastq \
#         --aligned ./sortmerna/aligned_barcode$i.porechopped \
#         --fastx --other ./sortmerna/unaligned_barcode$i.porechopped \
#         --blast '1 cigar qcov' barcode$i.blast.out
#         rm -rf /home/kebaso/sortmerna/run/
# done


# STEP 4: ISONCLUST & ISONCORRECT

# Should be run on unzipped fastq

#file=./sortmerna/*.porechopped.fq
#for i in ${file}
#do
#echo ${i}
 
#base=$(basename $i .porechopped.fq) 

#echo ${base}

#isONclust --ont --fastq $i --outfolder ./isONclust/$base.isonclust.output/

#isONclust write_fastq --clusters ./isONclust/$base.isonclust.output/final_clusters.tsv --fastq $i --outfolder ./isONclust/$base.isonclust.fastq.output --N 1

#move isoncorrect files to isonclust before running above.

#run_isoncorrect --fastq_folder ./isONclust/$base.isonclust.fastq.output --k 13 --w 20 --outfolder ./isONclust/$base.correction
#done

#OPTIONAL BELOW TO MERGE ALL CORRECTED READS INTO ONE FILE

# for outfolder in ./isONclust/*.correction; do 
#     #echo ${outfolder}
#     base=$(basename $outfolder .correction)
#     #echo ${base}
#     touch $outfolder/all_corrected_reads.fq

#     OUTFILES=$outfolder/*/corrected_reads.fastq
#     #echo ${OUTFILES}
    
#     for f in $OUTFILES; do 
#          echo $f
#          cat $f >> $outfolder/$base.all_corrected_reads.fq

# echo
# echo "Finished with pipeline and wrote corrected reads to: " $outfolder/$base.all_corrected_reads.fq
# echo
# done
# done

#copy corrected reads to one folder

# mkdir all_corrected_reads

# for outfolder in ./isONclust/*.correction; do 
# #echo $outfolder
#      base=$(basename $outfolder .correction)
#      #echo ${base}
#      cp $outfolder/$base.all_corrected_reads.fq all_corrected_reads/
#      echo "Finished copying corrected reads" $base.all_corrected_reads.fq "to: all_corrected_reads/"
# done


#Minmap2 
#Mapping to host genome and counting

#for i in all_corrected_reads/*.all_corrected_reads.fq; do 
#echo ${i}
 #      base=$(basename $i)
#       minimap2 -a -x splice -k13 AR.fa ${i} -o output_reads/${base}.aln.sam
#done

#for i in ./output_reads/*.aln.sam; do 
    #base=$(basename $i .aln.sam)
    #samtools view -bS $i > $base.bam
#done

#for i in *.bam; do 
#   base=$(basename $i .bam)
#   samtools view -h -F 4 $i -o $base.mapped.bam
#   samtools sort $base.mapped.bam > $base.mapped.sorted.bam
#   samtools index $base.mapped.sorted.bam
#done


#for i in *.mapped.sorted.bam; do 
#    base=$(basename $i .mapped.sorted.bam)
#    NanoCount -i $i -o $base.nanocount.tsv -b $base.nanocount.bam -t 0.8 --verbose
#done   

#Download reference transcriptome:
#wget https://www.ncbi.nlm.nih.gov/nuccore/NC_053517.1
#wget https://www.ncbi.nlm.nih.gov/nuccore/NC_053518.1
#wget https://www.ncbi.nlm.nih.gov/nuccore/NC_053519.1
#concatenate the files into one reference sequence.

#Mapping to host transcriptome and counting
#for i in ./all_corrected_reads/*all_corrected_reads.fq; do 
#    base=$(basename $i .all_corrected_reads.fq)
    #echo ${base}


#       minimap2 -a -x map-ont -k13 GCF_016920715.1_AaraD3_rna.fna $i -o output_reads_rna/$base.aln.rna.sam
#done