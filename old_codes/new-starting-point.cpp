#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

struct SamEntry {
    int flag;
    int pos;
    string cigar;
    string seq;
};

// Funkcija za parsiranje SAM datoteke
vector<SamEntry> parseSamFile(const string& filename) {
    vector<SamEntry> entries;
    ifstream infile(filename);
    if (!infile) {
        cerr << "Error opening file: " << filename << endl;
        return entries;
    }

    string line;
    while (getline(infile, line)) {
        if (line[0] == '@') continue; // preskoči zaglavlja

        istringstream iss(line);
        string qname, rname, cigar, rnext, seq, qual;
        int flag, pos, mapq, pnext, tlen;

        iss >> qname >> flag >> rname >> pos >> mapq >> cigar >> rnext >> pnext >> tlen >> seq >> qual;

        if (flag & 4) continue; // preskoči redove gde je FLAG 4

        SamEntry entry;
        entry.flag = flag;
        entry.pos = pos;
        entry.cigar = cigar;
        entry.seq = seq;

        entries.push_back(entry);
    }

    return entries;
}

// Funkcija za obradu CIGAR stringa
vector<pair<char, int>> parseCigar(const string& cigar) {
    vector<pair<char, int>> operations;
    istringstream iss(cigar);
    char op;
    int count;
    while (iss >> count >> op) {
        operations.push_back({op, count});
    }
    return operations;
}

// Funkcija za detekciju mutacija
void detectMutations(const vector<SamEntry>& samEntries, const string& reference, const string& outputCsv) {
    ofstream outfile(outputCsv);
    if (!outfile) {
        cerr << "Error opening file: " << outputCsv << endl;
        return;
    }

    outfile << "type,X,pos,base\n";

    map<int, map<char, int>> substitutionVotes;
    map<int, map<string, int>> insertionVotes;
    map<int, map<string, int>> deletionVotes;

    for (const auto& entry : samEntries) {
        int refPos = entry.pos - 1; // SAM format uses 1-based indexing
        int readPos = 0;
        string readSeq = entry.seq;

        // add value of first soft clipping (it is always the first element of cigar) to readPos
        int firstSoftClipping = 0;
        istringstream iss(entry.cigar);
        iss >> firstSoftClipping;

        cout << "Current soft clipping: " << firstSoftClipping << endl;

        
        // Ispis trenutne sekvence (prvih 50 znakova) 
        cout << "Read: ";
        for (int i = 0; i < 50; ++i) {
            if (readPos + i < readSeq.size()) {
                cout << readSeq[readPos + firstSoftClipping + i];
            }
        }

        cout << endl;
        
        vector<char> refOutput =  vector<char>(reference.begin() + refPos, reference.begin() + refPos + 50);
        
        cout << "Reference: ";
        for (char c : refOutput) {
            cout << c;
        }

        cout << endl;




        if (entry.flag & 16) {
            // Reverzni komplement, preokreni sekvencu
            reverse(readSeq.begin(), readSeq.end());
        }

        auto cigarOperations = parseCigar(entry.cigar);
        for (const auto& op_count : cigarOperations) {
            char op = op_count.first;
            int count = op_count.second;

            if (op == 'M') { // Match or mismatch
                for (int i = 0; i < count; ++i) {
                    if (refPos >= 0 && refPos < reference.size() && readPos >= 0 && readPos < readSeq.size() &&
                        reference[refPos] != readSeq[readPos]) {
                        substitutionVotes[refPos][readSeq[readPos]]++;
                    }
                    ++refPos;
                    ++readPos;
                }
            } else if (op == 'I') { // Insertion
                if (readPos >= 0 && readPos + count <= readSeq.size()) {
                    string insertedBases = readSeq.substr(readPos, count);
                    insertionVotes[refPos][insertedBases]++;
                    readPos += count;
                }
            } else if (op == 'D') { // Deletion
                if (refPos >= 0 && refPos + count <= reference.size()) {
                    string deletedBases = reference.substr(refPos, count);
                    deletionVotes[refPos][deletedBases]++;
                    refPos += count;
                }
            } else if (op == 'S') { // Soft clipping
                readPos += count;
            } else if (op == 'H') { // Hard clipping
                continue;
            }
        }
    }

    // Većinsko glasanje za supstitucije
    for (map<int, map<char, int>>::iterator it = substitutionVotes.begin(); it != substitutionVotes.end(); ++it) {
        int pos = it->first;
        map<char, int>& votes = it->second;

        map<char, int>::iterator maxVote = max_element(votes.begin(), votes.end(), [](const pair<char, int>& a, const pair<char, int>& b) {
            return a.second < b.second;
        });

        if (maxVote != votes.end() && maxVote->first != reference[pos]) {
            outfile << "substitution,X," << pos << "," << maxVote->first << "\n";
        }
    }

    // Većinsko glasanje za umetanja
    for (map<int, map<string, int>>::iterator it = insertionVotes.begin(); it != insertionVotes.end(); ++it) {
        int pos = it->first;
        map<string, int>& votes = it->second;

        map<string, int>::iterator maxVote = max_element(votes.begin(), votes.end(), [](const pair<string, int>& a, const pair<string, int>& b) {
            return a.second < b.second;
        });

        if (maxVote != votes.end()) {
            outfile << "insertion,I," << pos << "," << maxVote->first << "\n";
        }
    }

    // Većinsko glasanje za brisanja
    for (map<int, map<string, int>>::iterator it = deletionVotes.begin(); it != deletionVotes.end(); ++it) {
        int pos = it->first;
        map<string, int>& votes = it->second;

        map<string, int>::iterator maxVote = max_element(votes.begin(), votes.end(), [](const pair<string, int>& a, const pair<string, int>& b) {
            return a.second < b.second;
        });

        if (maxVote != votes.end()) {
            outfile << "deletion,D," << pos << "," << maxVote->first << "\n";
        }
    }

    outfile.close();
}

int main() {
    string samFile = "minimap_output.sam";
    string referenceFile = "lambda.fasta";
    string outputCsv = "mutations.csv";

    // Parsiraj SAM datoteku
    vector<SamEntry> samEntries = parseSamFile(samFile);

    // Učitaj referentni genom
    ifstream refFile(referenceFile);
    string reference;
    string line;
    while (getline(refFile, line)) {
        if (line[0] != '>') {
            reference += line;
        }
    }

    // Detekcija mutacija
    detectMutations(samEntries, reference, outputCsv);

    cout << "Mutation detection completed, output saved to " << outputCsv << endl;
    return 0;
}
