#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>

using namespace std;

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


int main() {
    // This is hard-coded for simplicity, but could be passed as command-line arguments
    string reference = "lambda.fasta";
    string reads = "lambda_simulated_reads.fasta";
    string samOutput = "alignment.sam";
    string csvOutput = "mutations.csv";

    // Step 1: Run Minimap2
    string minimapCmd = "minimap2 -ax map-ont " + reference + " " + reads + " > " + samOutput;
    cout << "Running Minimap2: " << minimapCmd << endl;
    executeCommand(minimapCmd);

    // Step 2: Parse SAM file to detect mutations
    detectMutations(samOutput, csvOutput);

    cout << "Mutation detection completed. Results are stored in " << csvOutput << endl;

    return 0;
}
