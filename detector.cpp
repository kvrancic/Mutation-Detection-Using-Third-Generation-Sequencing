#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;

// Structure to represent a mutation
struct Mutation {
    string Type;           // Type of mutation (Substitution, Insertion, Deletion, etc.)
    int Position;          // Position of the mutation
    string Base;           // Base or sequence involved in the mutation
};


// Function to execute system command
void executeCommand(const string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Command failed: " << command << endl;
        exit(EXIT_FAILURE);
    }
}

// Function to parse the SAM file and extract mutations
void detectMutations(const string& samFile, const string& outputFile) {
    ifstream inFile(samFile);
    ofstream outFile(outputFile);
    string line;

    // Write the header for the CSV file
    outFile << "Type,Position,Base\n";

    if (!inFile.is_open() || !outFile.is_open()) {
        cerr << "Failed to open files." << endl;
        exit(EXIT_FAILURE);
    }

    regex cigarRegex("([0-9]+)([MIDNSHP=X])");
    smatch matches;

    while (getline(inFile, line)) {
        if (line[0] == '@') continue; // Skip header lines

        istringstream iss(line);
        vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());

        if (tokens.size() > 9 && tokens[1] != "4") { // Aligned sequence
            string seq = tokens[9];
            int pos = stoi(tokens[3]); // Starting position of alignment
            string cigar = tokens[5];
            size_t seqIdx = 0;

            while (regex_search(cigar, matches, cigarRegex)) {
                int len = stoi(matches[1]);
                char type = matches[2].str()[0];

                switch (type) {
                    case 'M': // Match or mismatch (assumed match here)
                        seqIdx += len;
                        pos += len;
                        break;
                    case 'X': // Sequence mismatch explicitly
                        for (int i = 0; i < len; ++i) {
                            outFile << "Supstitucija,X," << (pos++) << "," << seq[seqIdx++] << "\n";
                        }
                        break;
                    case 'I': // Insertion to the reference
                        outFile << "Umetanje,I," << pos << ",";
                        for (int i = 0; i < len; ++i) {
                            outFile << seq[seqIdx++];
                        }
                        outFile << "\n";
                        break;
                    case 'D': // Deletion from the reference
                        outFile << "Brisanje,D," << pos << ",-\n";
                        pos += len;
                        break;
                    case 'N': // Skipped region from the reference
                        pos += len;
                        break;
                    case 'S': // Soft clipping
                    case 'H': // Hard clipping
                    case 'P': // Padding
                    case '=': // Sequence match
                        pos += len;
                        seqIdx += len;
                        break;
                }
                cigar = matches.suffix().str();
            }
        }
    }

    inFile.close();
    outFile.close();
}

// Function to load variants from the VCF file generated from FreeBayes
vector<Mutation> loadVariantsFromVCF(const string& vcfFile) {
    ifstream inFile(vcfFile);
    vector<Mutation> variants;
    string line;

    while (getline(inFile, line)) {
        if (line[0] == '#') continue; // Skip header lines

        istringstream iss(line);
        string chr, pos, id, ref, alt, qual, filter, info, format, sample;

        // Extract fields from the VCF line
        iss >> chr >> pos >> id >> ref >> alt >> qual >> filter >> info >> format >> sample;

        if (alt != ".") { // If variant already exists
            // Take only the first variant
            istringstream altStream(alt);
            string firstAlt;
            getline(altStream, firstAlt, ',');
            variants.push_back({"Variant", stoi(pos) - 1, firstAlt}); // -1 because of zero-based position
        }
    }

    return variants;
}


int main() {
    // This is hard-coded for simplicity, but could be passed as command-line arguments
    string reference = "data/lambda.fasta";
    string reads = "data/lambda_simulated_reads.fasta";
    string samOutput = "data/alignment.sam";
    string bamOutput = "data/alignment.bam";
    string sortedBam = "data/alignment_sorted.bam";
    string csvOutput = "data/mutations.csv";
    string vcfOutput = "data/variants.vcf";
    string indexCmd = "samtools faidx " + reference; // Command to generate index
    // Define paths for marked BAM file and metrics file
    string markedBam = "data/alignment_marked.bam"; 
    string metricsFile = "data/mark_duplicates_metrics.txt"; 


    // Step 1: Run Minimap2
    string minimapCmd = "minimap2 -ax map-ont " + reference + " " + reads + " > " + samOutput;
    cout << "Running Minimap2: " << minimapCmd << endl;
    executeCommand(minimapCmd);

    // Step 2: Parse SAM file to detect mutations
    detectMutations(samOutput, csvOutput);

    cout << "Mutation detection completed. Results are stored in " << csvOutput << endl;

    // Step 3: Generate index for the reference genome
    cout << "Generating index for the reference genome..." << endl;
    executeCommand(indexCmd); // Execute the command to generate index

    // Step 4: Convert SAM to BAM
    string samToBamCmd = "samtools view -bS " + samOutput + " > " + bamOutput;
    cout << "Converting SAM to BAM: " << samToBamCmd << endl;
    executeCommand(samToBamCmd);

    // Step 5: Sort BAM file
    string sortCmd = "samtools sort -o " + sortedBam + " " + bamOutput;
    cout << "Sorting BAM file: " << sortCmd << endl;
    executeCommand(sortCmd);

    // Step 6: Mark duplicates‚
    string markDuplicatesCmd = "picard MarkDuplicates I=" + sortedBam + " O=" + markedBam + " M=" + metricsFile + " REMOVE_DUPLICATES=true";
    cout << "Marking duplicates: " << markDuplicatesCmd << endl;
    executeCommand(markDuplicatesCmd);

    // Step 7: Index marked BAM file
    string indexMarkedCmd = "samtools index " + markedBam;
    cout << "Indexing marked BAM file: " << indexMarkedCmd << endl;
    executeCommand(indexMarkedCmd);


    // Step 8: Run FreeBayes for evaluation with optimization
    string freebayesCmd = "freebayes -f " + reference + " --use-best-n-alleles 4 --min-alternate-count 3 " + markedBam + " > " + vcfOutput;
    cout << "Running FreeBayes with optimization: " << freebayesCmd << endl;
    executeCommand(freebayesCmd);



    cout << "Variant calling completed. Results are stored in " << vcfOutput << endl;
     //inspect VCF file: bcftools view data/variants.vcf
 
    return 0;
}
