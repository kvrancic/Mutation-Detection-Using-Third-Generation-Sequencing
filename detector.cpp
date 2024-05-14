#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;

// Structure to store mutation information
struct Mutation {
    char type;
    int position;
    char nucleotide;

    Mutation(char t, int pos, char nuc) : type(t), position(pos), nucleotide(nuc) {}
};




// Function to execute system command
void executeCommand(const string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Command failed: " << command << endl;
        exit(EXIT_FAILURE);
    }
}


// Function to parse the CIGAR string and identify mutations
void parseCigarString(const string& cigar, int referenceStart, vector<Mutation>& mutations) {
    int position = referenceStart; // Initialize position on the reference genome

    // Parse the CIGAR string
    size_t pos = 0;
    while (pos < cigar.length()) {
        // Extract the length of each operation
        size_t end = cigar.find_first_of("MIDNSHP=X", pos);
        if (end == string::npos) {
            end = cigar.length();
        }
        int length = stoi(cigar.substr(pos, end - pos));

        // Determine the operation
        char operation = cigar[end];
        if (operation == 'X') {
            // Substitution mutation
            for (int i = 0; i < length; ++i) {
                int poz;
                char nucleotide; // Fetch the nucleotide from the read (for example)
                // Store the mutation information
                poz = position + i;
                mutations.push_back(Mutation('X', poz, nucleotide));
            }
        } else if (operation == 'I') {
            // Insertion mutation
            char nucleotide; // Fetch the inserted nucleotide from the read (for example)
            // Store the mutation information
            mutations.push_back(Mutation('I', position - 1, nucleotide)); // Position before insertion
        } else if (operation == 'D') {
            // Deletion mutation
            // Store the mutation information
            mutations.push_back(Mutation('D', position, '-'));
            // Move the reference position forward by the length of the deletion
            position += length;
        } else if (operation == 'M') {
            // Match operation, move the reference position forward
            position += length;
        }

        // Move to the next operation in the CIGAR string
        pos = end + 1;
    }
}

// Updated function to identify mutations from the PAF file
vector<Mutation> identifyMutationsFromPAF(const string& pafFilename) {
    vector<Mutation> mutations;

    ifstream file(pafFilename);
    if (!file.is_open()) {
        cerr << "Unable to open file " << pafFilename << endl;
        return mutations;
    }

    // Iterating through each line of the PAF file
    string line;
    while (getline(file, line)) {
        istringstream iss(line);
        string cigar;
        int referenceStart, referenceEnd;
        // Reading necessary data from the PAF line
        iss >> cigar >> referenceStart >> referenceEnd;
        // Parse the CIGAR string and identify mutations
        parseCigarString(cigar, referenceStart, mutations);
    }

    file.close();
    return mutations;
}


// Function to write mutations to a CSV file
void writeMutationsToCSV(const vector<Mutation>& mutations, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Unable to open file " << filename << endl;
        return;
    }

    // Writing CSV file header
    file << "Mutation,Position,Nucleotide" << endl;

    // Writing each mutation in CSV format
    for (const auto& mutation : mutations) {
        file << mutation.type << "," << mutation.position << "," << mutation.nucleotide << endl;
    }

    file.close();
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


int main() {
    // This is hard-coded for simplicity, but could be passed as command-line arguments
    string reference = "data/lambda.fasta";
    string reads = "data/lambda_simulated_reads.fasta";
    string samOutput = "data/alignment.sam";
    string pafOutput = "data/aligment2.paf";
    string bamOutput = "data/alignment.bam";
    string sortedBam = "data/alignment_sorted.bam";
    string csvOutput = "data/mutations.csv";
    string vcfOutput = "data/variants.vcf";
    string indexCmd = "samtools faidx " + reference; // Command to generate index
    // Define paths for marked BAM file and metrics file
    string markedBam = "data/alignment_marked.bam"; 
    string metricsFile = "data/mark_duplicates_metrics.txt"; 


    // Step 1: Run Minimap2
    string minimapCmd = "minimap2 --eqx -c -ax map-ont " + reference + " " + reads + " > " + pafOutput;
    cout << "Running Minimap2: " << minimapCmd << endl;
    executeCommand(minimapCmd);


    // Step 2: Identify mutations from the PAF file
    vector<Mutation> mutations = identifyMutationsFromPAF(pafOutput);
    // Step 3: Write mutations to a CSV file
    writeMutationsToCSV(mutations, csvOutput);

    cout << "Mutations are written to " << csvOutput << endl;

    // // Step 2: Parse SAM file to detect mutations
    // detectMutations(samOutput, csvOutput);

    // cout << "Mutation detection completed. Results are stored in " << csvOutput << endl;

    // // Step 3: Generate index for the reference genome
    // cout << "Generating index for the reference genome..." << endl;
    // executeCommand(indexCmd); // Execute the command to generate index

    // // Step 4: Convert SAM to BAM
    // string samToBamCmd = "samtools view -bS " + samOutput + " > " + bamOutput;
    // cout << "Converting SAM to BAM: " << samToBamCmd << endl;
    // executeCommand(samToBamCmd);

    // // Step 5: Sort BAM file
    // string sortCmd = "samtools sort -o " + sortedBam + " " + bamOutput;
    // cout << "Sorting BAM file: " << sortCmd << endl;
    // executeCommand(sortCmd);

    // // Step 6: Mark duplicates‚
    // string markDuplicatesCmd = "picard MarkDuplicates I=" + sortedBam + " O=" + markedBam + " M=" + metricsFile + " REMOVE_DUPLICATES=true";
    // cout << "Marking duplicates: " << markDuplicatesCmd << endl;
    // executeCommand(markDuplicatesCmd);

    // // Step 7: Index marked BAM file
    // string indexMarkedCmd = "samtools index " + markedBam;
    // cout << "Indexing marked BAM file: " << indexMarkedCmd << endl;
    // executeCommand(indexMarkedCmd);


    // // Step 8: Run FreeBayes for evaluation with optimization
    // string freebayesCmd = "freebayes -f " + reference + " --use-best-n-alleles 4 --min-alternate-count 3 " + markedBam + " > " + vcfOutput;
    // cout << "Running FreeBayes with optimization: " << freebayesCmd << endl;
    // //executeCommand(freebayesCmd);



    // cout << "Variant calling completed. Results are stored in " << vcfOutput << endl;
    //  //inspect VCF file: bcftools view data/variants.vcf
 
    return 0;
}
