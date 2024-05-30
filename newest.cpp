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

    struct MutationProposal {
        int substitutionVotes;
        char substitutionBase;
        int insertionVotes;
        string insertionBases;
        int deletionVotes;
        string deletionBases;
    };

    map<int, MutationProposal> mutationProposals;

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
                    if (refPos >= 0 && refPos < reference.size() && readPos >= 0 && readPos < readSeq.size() &&
                        reference[refPos] != readSeq[readPos]) {
                        mutationProposals[refPos].substitutionVotes++;
                        mutationProposals[refPos].substitutionBase = readSeq[readPos];
                    }
                    ++refPos;
                    ++readPos;
                }
            } else if (op == 'I') { // Insertion
                if (readPos >= 0 && readPos + count <= readSeq.size()) {
                    string insertedBases = readSeq.substr(readPos, count);
                    mutationProposals[refPos].insertionVotes++;
                    mutationProposals[refPos].insertionBases = insertedBases;
                    readPos += count;
                }
            } else if (op == 'D') { // Deletion
                if (refPos >= 0 && refPos + count <= reference.size()) {
                    string deletedBases = reference.substr(refPos, count);
                    mutationProposals[refPos].deletionVotes++;
                    mutationProposals[refPos].deletionBases = deletedBases;
                    refPos += count;
                }
            } else if (op == 'S') { // Soft clipping
                readPos += count;
            } else if (op == 'H') { // Hard clipping
                continue;
            }
        }
    }

    // Većinsko glasanje za svaku poziciju
    for (auto& proposal : mutationProposals) {
        int pos = proposal.first;
        MutationProposal& votes = proposal.second;

        int maxVotes = max({votes.substitutionVotes, votes.insertionVotes, votes.deletionVotes});
        if (maxVotes == votes.substitutionVotes && votes.substitutionBase != reference[pos]) {
            outfile << "substitution,X," << pos << "," << votes.substitutionBase << "\n";
        } else if (maxVotes == votes.insertionVotes) {
            outfile << "insertion,I," << pos << "," << votes.insertionBases << "\n";
        } else if (maxVotes == votes.deletionVotes) {
            outfile << "deletion,D," << pos << "," << votes.deletionBases << "\n";
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
