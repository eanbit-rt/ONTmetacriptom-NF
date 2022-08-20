#!/usr/bin/env bash

# Download the sortmerna package into rRNA_databases folder
wget --quiet \
     --directory-prefix=rRNA_databases \
     https://github.com/biocore/sortmerna/archive/2.1b.zip

# Decompress folder 
unzip -q rRNA_databases/2.1b.zip -d rRNA_databases

# Remove the zipped file
rm -rf rRNA_databases/*.zip

# Move the database files into the correct directory
mv rRNA_databases/sortmerna-2.1b/rRNA_databases/silva*.fasta rRNA_databases/

# Remove the Ddecompressed directory
rm -rf rRNA_databases/sortmerna*
