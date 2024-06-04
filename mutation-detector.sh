#!/bin/bash

set -e

# Hard-coded paths for simplicity
reference="data/lambda.fasta"
reads="data/lambda_simulated_reads.fasta"
samOutput="data/alignment.sam"
bamOutput="data/alignment.bam"
sortedBam="data/alignment_sorted.bam"
csvOutput="data/mutations.csv"
vcfOutput="data/variants.vcf"
markedBam="data/alignment_marked.bam"

# Step 1: Run Minimap2
minimapCmd="minimap2 -ax map-ont $reference $reads > $samOutput"
echo "Running Minimap2: $minimapCmd"
eval $minimapCmd

# Step 2: Call detector to detect mutations
# it is a C++ program that takes the SAM file and reference genome as input, so please first convert .cpp to exe 
# write that as a part of the script 

detectorCmd="./detector $samOutput $reference $csvOutput"
echo "Running mutation detection: $detectorCmd"
eval $detectorCmd

echo "Mutation detection completed. Results are stored in $csvOutput"

# Step 3: Generate index for the reference genome
indexCmd="samtools faidx $reference"
echo "Generating index for the reference genome..."
eval $indexCmd

# Step 4: Convert SAM to BAM
samToBamCmd="samtools view -bS $samOutput > $bamOutput"
echo "Converting SAM to BAM: $samToBamCmd"
eval $samToBamCmd

# Step 5: Sort BAM file
sortCmd="samtools sort -o $sortedBam $bamOutput"
echo "Sorting BAM file: $sortCmd"
eval $sortCmd

# Step 6: Index marked BAM file
indexMarkedCmd="samtools index $markedBam"
echo "Indexing marked BAM file: $indexMarkedCmd"
eval $indexMarkedCmd

# Step 7: Run FreeBayes for evaluation with optimization
freebayesCmd="freebayes -f $reference --use-best-n-alleles 4 --min-alternate-count 3 $markedBam > $vcfOutput"
echo "Running FreeBayes with optimization: $freebayesCmd"
eval $freebayesCmd

echo "Variant calling completed. Results are stored in $vcfOutput"
