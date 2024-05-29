#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

// Struktura za pohranu informacija iz SAM datoteke
struct SamEntry {
    string qname;
    int flag;
    string rname;
    int pos;
    int mapq;
    string cigar;
    string seq;
    string qual;
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

        SamEntry entry;
        istringstream iss(line);
        iss >> entry.qname >> entry.flag >> entry.rname >> entry.pos >> entry.mapq >> entry.cigar >> entry.seq >> entry.qual;

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

    outfile << "type,X,pos,reference,read\n";

    int windowSize = 1; // širina prozora za provjeru poravnanja

    for (const auto& entry : samEntries) {
        int refPos = entry.pos - 1; // SAM format uses 1-based indexing
        int readPos = 0;
        string readSeq = entry.seq;

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
                    if (refPos >= reference.size() || readPos >= readSeq.size()) {
                        break;
                    }

                    if (reference[refPos] == readSeq[readPos]) {
                        // Ako se poklapaju, idemo dalje
                        ++refPos;
                        ++readPos;
                    } else {
                        // Ako se ne poklapaju, proširujemo prozor i provjeravamo
                        bool found = false;
                        for (int w = 1; w <= windowSize; ++w) {
                            // Provjera umetanja
                            if (readPos + w < readSeq.size() && reference[refPos] == readSeq[readPos + w]) {
                                outfile << "insertion,I," << refPos << "," << readSeq.substr(readPos, w) << "\n";
                                readPos += w + 1; // Preskačemo umetnute elemente
                                ++refPos;
                                found = true;
                                break;
                            }
                            // Provjera brisanja
                            if (refPos + w < reference.size() && reference[refPos + w] == readSeq[readPos]) {
                                outfile << "deletion,D," << refPos << "," << reference.substr(refPos, w) << "\n";
                                refPos += w + 1; // Preskačemo obrisane elemente
                                ++readPos;
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            // Ako se ni jedna mutacija ne poklapa, bilježimo supstituciju
                            outfile << "substitution,X," << refPos << "," << reference[refPos] << "," << readSeq[readPos] << "\n";
                            ++refPos;
                            ++readPos;
                        }
                    }
                }
            } else if (op == 'I') { // Insertion
                if (readPos >= 0 && readPos + count <= readSeq.size()) {
                    string insertedBases = readSeq.substr(readPos, count);
                    outfile << "insertion,I," << refPos << "," << insertedBases << "\n";
                    readPos += count;
                }
            } else if (op == 'D') { // Deletion
                if (refPos >= 0 && refPos + count <= reference.size()) {
                    string deletedBases = reference.substr(refPos, count);
                    outfile << "deletion,D," << refPos << "," << deletedBases << "\n";
                    refPos += count;
                }
            } else if (op == 'S') { // Soft clipping
                readPos += count;
            } else if (op == 'H') { // Hard clipping
                continue;
            }
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
