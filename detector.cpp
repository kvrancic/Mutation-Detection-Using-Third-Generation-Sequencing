#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <regex>
#include <algorithm>

using namespace std;

struct Mutation {
    string type;  // Type of mutation (Substitution, Insertion, Deletion)
    int position; // Position of the mutation
    string base;  // Base or sequence involved in the mutation
};

void executeCommand(const string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Command failed: " << command << endl;
        exit(EXIT_FAILURE);
    }
}

void mapReadsWithMinimap2(const string& reference, const string& reads, const string& samOutput) {
    string minimapCmd = "minimap2 -ax map-ont " + reference + " " + reads + " > " + samOutput;
    cout << "Running Minimap2: " << minimapCmd << endl;
    executeCommand(minimapCmd);
}

void detectMutations(const string& samFile, const string& referenceFile, const string& outputFile) {
    ifstream sam(samFile);
    ifstream reference(referenceFile);
    ofstream outFile(outputFile);
    string line;

    // Read reference genome into a string
    string refGenome, refLine;
    while (getline(reference, refLine)) {
        if (refLine[0] == '>') continue; // Skip the header line
        refGenome += refLine;
    }

    // Write the header for the CSV file
    outFile << "Type,Position,Base\n";

    if (!sam.is_open() || !outFile.is_open()) {
        cerr << "Failed to open files." << endl;
        exit(EXIT_FAILURE);
    }

    map<int, map<char, int>> positionBaseCount;
    regex cigarRegex("([0-9]+)([MIDNSHP=X])");
    smatch matches;

    while (getline(sam, line)) {
        if (line[0] == '@') continue; // Skip header lines

        istringstream iss(line);
        vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());

        if (tokens.size() > 9 && tokens[1] != "4") { // Aligned sequence
            string seq = tokens[9];
            int pos = stoi(tokens[3]) - 1; // Starting position of alignment (0-based index)
            string cigar = tokens[5];
            size_t seqIdx = 0;

            while (regex_search(cigar, matches, cigarRegex)) {
                int len = stoi(matches[1]);
                char type = matches[2].str()[0];

                switch (type) {
                    case 'M': // Match or mismatch
                        for (int i = 0; i < len; ++i) {
                            char refBase = refGenome[pos + i];
                            char readBase = seq[seqIdx + i];
                            if (refBase != readBase) {
                                positionBaseCount[pos + i][readBase]++;
                            }
                        }
                        pos += len;
                        seqIdx += len;
                        break;
                    case 'I': // Insertion
                        for (int i = 0; i < len; ++i) {
                            positionBaseCount[pos]['I']++;
                        }
                        seqIdx += len;
                        break;
                    case 'D': // Deletion
                        for (int i = 0; i < len; ++i) {
                            positionBaseCount[pos + i]['D']++;
                        }
                        pos += len;
                        break;
                    case 'S': // Soft clipping
                    case 'H': // Hard clipping
                    case 'N': // Skipped region from the reference
                    case 'P': // Padding
                    case '=': // Sequence match
                    case 'X': // Sequence mismatch
                        break;
                }
                cigar = matches.suffix().str();
            }
        }
    }

    // Detect mutations based on majority vote
    for (const auto& entry : positionBaseCount) {
        int pos = entry.first;
        const auto& baseCounts = entry.second;

        if (baseCounts.size() > 1) {
            auto maxElement = max_element(baseCounts.begin(), baseCounts.end(),
                                          [](const pair<char, int>& a, const pair<char, int>& b) {
                                              return a.second < b.second;
                                          });

            char majorityBase = maxElement->first;
            if (majorityBase == 'I') {
                string insertedBases;
                for (const auto& baseCount : baseCounts) {
                    if (baseCount.first != 'I') {
                        insertedBases += baseCount.first;
                    }
                }
                outFile << "I," << pos << "," << insertedBases << "\n";
            } else if (majorityBase == 'D') {
                outFile << "D," << pos << ",-\n";
            } else {
                char refBase = refGenome[pos];
                if (refBase != majorityBase) {
                    outFile << "X," << pos << "," << majorityBase << "\n";
                }
            }
        }
    }

    sam.close();
    reference.close();
    outFile.close();
}

void convertSamToBam(const string& samFile, const string& bamFile) {
    string samToBamCmd = "samtools view -bS " + samFile + " > " + bamFile;
    cout << "Converting SAM to BAM: " << samToBamCmd << endl;
    executeCommand(samToBamCmd);
}

void sortBamFile(const string& bamFile, const string& sortedBamFile) {
    string sortCmd = "samtools sort -o " + sortedBamFile + " " + bamFile;
    cout << "Sorting BAM file: " << sortCmd << endl;
    executeCommand(sortCmd);
}

void indexBamFile(const string& bamFile) {
    string indexCmd = "samtools index " + bamFile;
    cout << "Indexing BAM file: " << indexCmd << endl;
    executeCommand(indexCmd);
}

void runFreeBayes(const string& reference, const string& bamFile, const string& vcfOutput) {
    string freebayesCmd = "freebayes -f " + reference + " " + bamFile + " > " + vcfOutput;
    cout << "Running FreeBayes: " << freebayesCmd << endl;
    executeCommand(freebayesCmd);
}

int main() {
    string reference = "data/lambda.fasta";
    string reads = "data/lambda_mutated.fasta";
    string samOutput = "data/alignment.sam";
    string bamOutput = "data/alignment.bam";
    string sortedBam = "data/alignment_sorted.bam";
    string csvOutput = "data/mutations.csv";
    string vcfOutput = "data/variants.vcf";

    // Step 1: Run Minimap2
    mapReadsWithMinimap2(reference, reads, samOutput);

    cout << "Mapping completed. SAM file is stored in " << samOutput << endl;

    // Step 2: Detect mutations
    detectMutations(samOutput, reference, csvOutput);

    cout << "Mutation detection completed. Results are stored in " << csvOutput << endl;

    // Step 3: Convert SAM to BAM
    convertSamToBam(samOutput, bamOutput);

    // Step 4: Sort BAM file
    sortBamFile(bamOutput, sortedBam);

    // Step 5: Index BAM file
    indexBamFile(sortedBam);

    // Step 6: Run FreeBayes
    runFreeBayes(reference, sortedBam, vcfOutput);

    cout << "Variant calling completed. Results are stored in " << vcfOutput << endl;

    return 0;
}
