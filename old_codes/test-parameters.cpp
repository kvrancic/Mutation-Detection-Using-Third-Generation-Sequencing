#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <iterator>
#include <algorithm>

using namespace std;

// Function to execute system command
void executeCommand(const string& command) {
    int result = system(command.c_str());
    if (result != 0) {
        cerr << "Command failed: " << command << endl;
        exit(EXIT_FAILURE);
    }
}

// Function to parse the SAM file and calculate the success rate
double calculateSuccessRate(const string& samFile) {
    ifstream inFile(samFile);
    string line;
    int totalSequences = 0;
    int alignedSequences = 0;

    if (!inFile.is_open()) {
        cerr << "Failed to open SAM file." << endl;
        exit(EXIT_FAILURE);
    }

    while (getline(inFile, line)) {
        if (line[0] == '@') continue; // Skip header lines

        totalSequences++;
        istringstream iss(line);
        vector<string> tokens((istream_iterator<string>(iss)), istream_iterator<string>());

        if (tokens.size() > 1 && tokens[2] != "*") { // Aligned sequence
            alignedSequences++;
        }
    }

    inFile.close();
    return static_cast<double>(alignedSequences) / totalSequences;
}

int main() {
    // This is hard-coded for simplicity, but could be passed as command-line arguments
    string reference = "data/lambda.fasta";
    string reads = "data/lambda_simulated_reads.fasta";
    string samOutput = "data/alignment.sam";
    string resultFile = "data/success_rates.txt"; // Path to the result file

    // Open the result file
    ofstream outFile(resultFile);
    if (!outFile.is_open()) {
        cerr << "Failed to open result file." << endl;
        return EXIT_FAILURE;
    }

    // Define different parameter sets for Minimap2
    vector<string> minimapParams = {
        "-ax map-ont",
        "-ax sr",
        "-ax splice",
        "-ax splice:hq",
        "-ax map-pb",
        "-ax asm5",
        "-ax asm10",
        "-ax ava-ont",
        "-ax ava-pb"
    };

    vector<double> successRates;

    for (const auto& params : minimapParams) {
        string minimapCmd = "minimap2 " + params + " " + reference + " " + reads + " > " + samOutput;
        outFile << "Running Minimap2 with params: " << params << endl;
        executeCommand(minimapCmd);

        double successRate = calculateSuccessRate(samOutput);
        successRates.push_back(successRate);

        outFile << "Success rate for params (" << params << "): " << successRate << endl;
    }

    // Find the best parameter set
    auto maxElementIter = max_element(successRates.begin(), successRates.end());
    int bestIndex = distance(successRates.begin(), maxElementIter);
    outFile << "Best parameter set: " << minimapParams[bestIndex] << " with success rate: " << *maxElementIter << endl;

    // Close the result file
    outFile.close();

    return 0;
}
