#!/usr/bin/env bash

path="${1}"
files="barcode*"
dirpath=${path}/${files}

#---------------- STEP 1: RUNNING nanoQC ----------------------------
# Get all directories stating with 'barcode'
for path_per_dir in ${dirpath}; do 
    dirname=$(basename ${path_per_dir})
    # Get all file per directory
    echo "STEP 1: Running nanoQC for ${dirname}..."
    for file_with_path in ${path_per_dir}/*; do
        filename=$(basename --suffix=.fastq.gz $file_with_path)
        nanoQC -o nanoQC_Ouput/$dirname ${file_with_path} &>/dev/null
        # Rename the ouput file
        mv nanoQC_Ouput/${dirname}/nanoQC.html nanoQC_Ouput/${dirname}/${filename}.html
    done
    echo "done!"
done
echo "STEP 1: Completed successfully!"

# -------------------------- end ----------------------------


# ----------- RUNNING porechop ------------------------------
# porechop removes adapters
# Combo script that runs porechop on all fastq files and 
# concatenates them into one porechopped fastq 
# Get all directories stating with 'barcode'
echo "STEP 2: Running porechop..."
for path_per_dir in ${dirpath}; do 
    dirname=$(basename ${path_per_dir})
    mkdir -p porechop_temp/${dirname} porechop_output
    # Get all file per directory
    for file_with_path in ${path_per_dir}/*; do
        filename=$(basename --suffix=.fastq.gz $file_with_path)
        porechop \
            -i  ${file_with_path} \
            --format fastq \
            -o porechop_temp/${dirname}/${filename}_porechopped.fastq
        # Concatinating
        cat \
            porechop_temp/${dirname}/${filename}_porechopped.fastq >> \
            porechop_output/${dirname}_porechopped_concatenated.fastq

        # remove the porechopped sequence 
        rm porechop_temp/${dirname}/${filename}_porechopped.fastq 
    done
done
rm -r porechop_temp
echo "STEP 2: Completed successfully!"
echo "Moved porechop results to ${PWD}/porechop_output"
# -------------------end ----------------------------------------------------------


# #------------------------- Downloading the rNA Databases -----------------
# First check wether the rRNA Databases directory exits
if [ -d rRNA_databases ] && [ ! -z "$(ls -A -- rRNA_databases)" ]; then
    echo "rRNA Databases availables in ${PWD}/rRNA_databases directory"
else
    # # Download the sortmerna package into rRNA_databases folder
    echo "Downloading rRNA Databases..."
    wget --quiet \
        --directory-prefix=rRNA_databases \
        https://github.com/biocore/sortmerna/archive/2.1b.zip

    # Decompress folder 
    echo "Unzipping rRNA Databases..."
    unzip -q rRNA_databases/2.1b.zip -d rRNA_databases

    echo "Tidying up things..."
    # Remove the zipped file
    rm -rf rRNA_databases/*.zip

    # Move the database files into the correct directory
    mv rRNA_databases/sortmerna-2.1b/rRNA_databases/silva*.fasta rRNA_databases/

    # Remove the Decompressed directory
    rm -rf rRNA_databases/sortmerna*

    echo "Done"
    echo "rRNA Databases availables in ${PWD}/rRNA_databases directory"
fi

# ----------------------------------end-------------------------------------


# ------------------------------- RUNNING sortmerna tool -----------------
# Sorts rRNA and mRNA using reference rRNA databases
porechopped_files="./porechop_output/*"
mkdir -p rRNABlastOutput rRNALogOutput
for file in ${porechopped_files}; do
    filename=$(basename --suffix=.fastq ${file})
    sortmerna \
    --ref ./rRNA_databases/silva-arc-16s-id95.fasta \
    --ref ./rRNA_databases/silva-arc-23s-id98.fasta \
    --ref ./rRNA_databases/silva-bac-16s-id90.fasta \
    --ref ./rRNA_databases/silva-bac-23s-id98.fasta \
    --ref ./rRNA_databases/silva-euk-18s-id95.fasta \
    --ref ./rRNA_databases/silva-euk-28s-id98.fasta \
    -reads ${file} \
    --aligned rRNA_sequences/${filename}_sorted_rRNA \
    --fastx --other mRNA_sequences/${filename}_sorted_mRNA \
    --blast '1 cigar qcov' ${filename} 

    mv ${PWD}/rRNA_sequences/*.blast ${PWD}/rRNABlastOutput
    mv ${PWD}/rRNA_sequences/*.log ${PWD}/rRNALogOutput

    rm -r ${HOME}/sortmerna/run/kvdb
done
echo "STEP 3: Completed successfully!"

# -------------------------------end ------------------------