#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cstdlib>

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
// Function to detect mutations from .SAM file
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

    while (getline(inFile, line)) {
        if (line[0] == '@') continue; // Skip header lines

        istringstream iss(line);
        vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());

        // Assuming fields are properly formatted and the sequence is in the 10th field (index 9)
        if (tokens.size() > 9 && tokens[1] != "4") { // Make sure it's an aligned sequence
            string seq = tokens[9];
            int position = stoi(tokens[3]); // Starting position of alignment

            // Placeholder for mutation detection logic
            for (size_t i = 0; i < seq.length(); i++) {
                // This is a simplification; real mutation detection should consider the reference sequence
                outFile << "Substitution," << (position + i) << "," << seq[i] << "\n";
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
    executeCommand(minimapCmd);

    // Step 2: Parse SAM file to detect mutations
    detectMutations(samOutput, csvOutput);

    cout << "Mutation detection completed. Results are stored in " << csvOutput << endl;

    return 0;
}
