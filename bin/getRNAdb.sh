#!/usr/bin/env bash

# Downloading the rNA Databases
# First check wether the rRNA Databases directory exits
# and it is not empty.
if [ -d ${PWD}/rRNA_databases ] && [ ! -z "$(ls -A -- ${PWD}/rRNA_databases)" ]; then
    echo "rRNA Databases availables in ${PWD}/rRNA_databases directory"
    exit 0
else
    # Otherwise, download the sortmerna package, clean it up
    # into rRNA_databases folder 
    echo "Downloading rRNA Databases..."
    wget \
        --quiet \
        --directory-prefix=rRNA_databases \
        https://github.com/biocore/sortmerna/archive/2.1b.zip
    
    if [ "${$?}" -eq 0 ]; then
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
        exit 0
    else
        echo "ERROR: Unable to download rRNAdatabase! Check your internet connection"
        exit 1
    fi
fi