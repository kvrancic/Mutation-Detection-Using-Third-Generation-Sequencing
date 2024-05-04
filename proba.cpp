#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

// Function for parsing SAM records and writing mutations to a CSV file
void processSAM(const string& samFileName, const string& csvFileName) {
    ifstream samFile(samFileName);
    ofstream csvFile(csvFileName);

    // Check if files are opened successfully
    if (!samFile.is_open() || !csvFile.is_open()) {
        cerr << "Failed to open files!" << endl;
        return;
    }

    // Iterate through each line in the SAM file
    string line;
    while (getline(samFile, line)) {
        // Skip header lines in the SAM file
        if (line[0] == '@') continue;

        // Split line into fields
        vector<string> fields;
        string field;
        stringstream ss(line);
        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // Extract necessary information from the SAM record
        string cigar = fields[5];           // CIGAR string
        string referencePos = fields[3];    // Reference position
        string readSeq = fields[9];         // Read sequence

        // Iterate through the CIGAR string and identify mutations
        size_t pos = 0;
        size_t refPos = stoi(referencePos);
        for (char c : cigar) {
            if (isdigit(c)) continue;
            char op = c;
            int length = stoi(cigar.substr(pos, string::npos));
            pos += to_string(length).size();

            // If a mutation is detected, write it to the CSV file
            if (op == 'D' || op == 'I' || op == 'X') {
                if (op == 'D') {
                    csvFile << "D;" << refPos << ";" << "-" << endl;
                } else if (op == 'I') {
                    csvFile << "I;" << refPos - 1 << ";" << readSeq.substr(0, length) << endl;
                } else if (op == 'X') {
                    csvFile << "X;" << refPos << ";" << readSeq.substr(0, 1) << ">" << readSeq.substr(1, 1) << endl;
                }
            }

            // Update reference position
            if (op != 'I') refPos += length;
        }
    }

    // Close files
    samFile.close();
    csvFile.close();
}

int main() {
    // Mapping reads to reference genomes
    //system("minimap2 -c ecoli.fasta ecoli_simulated_reads.fasta > ecoli_alignments.sam");
    //system("minimap2 -c lambda.fasta lambda_simulated_reads.fasta > lambda_alignments.sam");

    // Mutation analysis
    processSAM("data/ecoli_alignments.sam", "data/ecoli_mutated.csv");

    //processSAM("data/ecoli_alignments.sam", "data/ecoli_mutaded.csv");
   // processSAM("lambda_alignments.sam", "lambda_mutaded.csv");

    // Evaluation with FreeBayes (assumed steps)
    // system("freebayes -f ecoli.fasta ecoli_simulated_reads.fasta > ecoli_variants.vcf");
    // system("freebayes -f lambda.fasta lambda_simulated_reads.fasta > lambda_variants.vcf");

    cout << "CSV files with mutations have been generated successfully." << endl;
    return 0;
}
