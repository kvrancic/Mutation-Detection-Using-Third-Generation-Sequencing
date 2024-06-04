#!/bin/bash

set -e

# Parse command line arguments
if [ $# -eq 0 ]; then
    echo "Usage: $0 <genome>"
    echo "Available genomes: lambda, ecoli"
    exit 1
fi

genome="$1"

case $genome in
    lambda)
        reference="data/lambda.fasta"
        ref_results_file="data/lambda_mutated.csv"
        ;;
    ecoli)
        reference="data/ecoli.fasta"
        ref_results_file="data/ecoli_mutated.csv"
        ;;
    *)
        echo "Invalid genome. Please choose 'lambda' or 'ecoli'."
        exit 1
        ;;
esac

# Paths for other files
reads="data/${genome}_simulated_reads.fasta"
samOutput="data/${genome}_alignment.sam"
bamOutput="data/alignment.bam"
sortedBam="data/alignment_sorted.bam"
csvOutput="data/mutations.csv"
vcfOutput="data/variants.vcf"
markedBam="data/alignment_marked.bam"
log_file="evaluation_log.txt"

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
#eval $freebayesCmd

echo "Variant calling completed. Results are stored in $vcfOutput"

# Call evaluate.py with appropriate parameters
python3 evaluate.py $csvOutput $ref_results_file $log_file
