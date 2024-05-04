#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cstdlib>

using namespace std;

// Funkcija za parsiranje SAM zapisa i zapisivanje mutacija u CSV formatu
void processSAM(const string& samFileName, const string& csvFileName) {
    ifstream samFile(samFileName);
    ofstream csvFile(csvFileName);

    if (!samFile.is_open() || !csvFile.is_open()) {
        cerr << "Failed to open files!" << endl;
        return;
    }

    // Prolazak kroz svaki redak u SAM datoteci
    string line;
    while (getline(samFile, line)) {
        // Preskačemo zaglavlje SAM datoteke
        if (line[0] == '@') continue;

        // Razdvajanje redaka na polja
        vector<string> fields;
        string field;
        stringstream ss(line);
        while (getline(ss, field, '\t')) {
            fields.push_back(field);
        }

        // Uzimamo potrebne informacije iz SAM zapisa
        string cigar = fields[5];
        string referencePos = fields[3];
        string readSeq = fields[9];

        // Iteriramo kroz CIGAR string i identificiramo mutacije
        size_t pos = 0;
        size_t refPos = stoi(referencePos);
        for (char c : cigar) {
            if (isdigit(c)) continue;
            char op = c;
            int length = stoi(cigar.substr(pos, string::npos));
            pos += to_string(length).size();

            // Ako je detektirana mutacija, zapisujemo je u CSV datoteku
            if (op == 'D' || op == 'I' || op == 'X') {
                if (op == 'D') {
                    csvFile << "D;" << refPos << ";" << "-" << endl;
                } else if (op == 'I') {
                    csvFile << "I;" << refPos - 1 << ";" << readSeq.substr(0, length) << endl;
                } else if (op == 'X') {
                    csvFile << "X;" << refPos << ";" << readSeq.substr(0, 1) << ">" << readSeq.substr(1, 1) << endl;
                }
            }

            // Ažuriramo referentnu poziciju
            if (op != 'I') refPos += length;
        }
    }

    // Zatvaramo datoteke
    samFile.close();
    csvFile.close();
}

int main() {
    // Mapiranje očitanja na referentne genome
    //system("minimap2 -c ecoli.fasta ecoli_simulated_reads.fasta > ecoli_alignments.sam");
    //system("minimap2 -c lambda.fasta lambda_simulated_reads.fasta > lambda_alignments.sam");

    // Analiza mutacija
    processSAM("data/ecoli_alignments.sam", "data/ecoli_mutated.csv");

    //processSAM("data/ecoli_alignments.sam", "data/ecoli_mutaded.csv");
   // processSAM("lambda_alignments.sam", "lambda_mutaded.csv");

    // Evaluacija s FreeBayesom (pretpostavljeni koraci)
    // system("freebayes -f ecoli.fasta ecoli_simulated_reads.fasta > ecoli_variants.vcf");
    // system("freebayes -f lambda.fasta lambda_simulated_reads.fasta > lambda_variants.vcf");

    cout << "CSV files with mutations have been generated successfully." << endl;
    return 0;
}
