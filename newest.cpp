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

// Funkcija za ispis vektora char-ova
ostream& operator<<(ostream& os, const vector<char>& vec) {
    for (char c : vec) {
        os << c;
    }
    return os;
}

// Funkcija za ispis vektora stringova
// ostream& operator<<(ostream& os, const vector<string>& vec) {
//     for (const string& str : vec) {
//         os << str << " ";
//     }
//     return os;
// }


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
        operations.push_back(std::make_pair(op, count));
    }
    return operations;
}

struct MutationProposal {
    int substitutionVotes = 0;
    vector<char> substitutionBases;
    int insertionVotes = 0;
    vector<string> insertionBases;
    int deletionVotes = 0;
    vector<string> deletionBases;
};

// fja za dodavanje na sve pozicije a ne samo pocetnu
void addMutationProposal(map<int, MutationProposal>& mutationProposals, int startPos, int count, char mutationType, const string& bases) {
    for (int i = 0; i < count; ++i) {
        int pos = startPos + i;
        if (mutationType == 'M') {
            mutationProposals[pos].substitutionVotes++;
            mutationProposals[pos].substitutionBases.push_back(bases[i]);
        } else if (mutationType == 'I') {
            mutationProposals[pos].insertionVotes++;
            mutationProposals[pos].insertionBases.push_back(bases.substr(i, 1));
        } else if (mutationType == 'D') {
            mutationProposals[pos].deletionVotes++;
            mutationProposals[pos].deletionBases.push_back(bases.substr(i, 1));
        }
    }
}

// Funkcija za detekciju mutacija
void detectMutations(const vector<SamEntry>& samEntries, const string& reference, const string& outputCsv) {
    ofstream outfile(outputCsv);
    if (!outfile) {
        cerr << "Error opening file: " << outputCsv << endl;
        return;
    }

    outfile << "type,X,pos,base\n";

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
                        addMutationProposal(mutationProposals, refPos, 1, 'M', string(1, readSeq[readPos]));
                    }
                    ++refPos;
                    ++readPos;
                }
            } else if (op == 'I') { // Insertion
                if (readPos >= 0 && readPos + count <= readSeq.size()) {
                    string insertedBases = readSeq.substr(readPos, count);
                    addMutationProposal(mutationProposals, refPos, 1, 'I', insertedBases);
                    readPos += count;
                }
            } else if (op == 'D') { // Deletion
                if (refPos >= 0 && refPos + count <= reference.size()) {
                    string deletedBases = reference.substr(refPos, count);
                    addMutationProposal(mutationProposals, refPos, count, 'D', deletedBases);
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

        // Pronalaženje najčešćeg znaka za svaku vrstu mutacije
        char mostCommonSubstitutionBase = ' ';
        if (!votes.substitutionBases.empty()) {
            map<char, int> substitutionBaseCounts;
            for (char base : votes.substitutionBases) {
                substitutionBaseCounts[base]++;
            }
            mostCommonSubstitutionBase = max_element(substitutionBaseCounts.begin(), substitutionBaseCounts.end(),
                                                     [](const pair<char, int>& a, const pair<char, int>& b) {
                                                         return a.second < b.second;
                                                     })->first;
        }

        string mostCommonInsertionBase;
        if (!votes.insertionBases.empty()) {
            map<string, int> insertionBaseCounts;
            for (const string& bases : votes.insertionBases) {
                insertionBaseCounts[bases]++;
            }
            mostCommonInsertionBase = max_element(insertionBaseCounts.begin(), insertionBaseCounts.end(),
                                                   [](const pair<string, int>& a, const pair<string, int>& b) {
                                                       return a.second < b.second;
                                                   })->first;
        }

        string mostCommonDeletionBase;
        if (!votes.deletionBases.empty()) {
            map<string, int> deletionBaseCounts;
            for (const string& bases : votes.deletionBases) {
                deletionBaseCounts[bases]++;
            }
            mostCommonDeletionBase = max_element(deletionBaseCounts.begin(), deletionBaseCounts.end(),
                                                 [](const pair<string, int>& a, const pair<string, int>& b) {
                                                     return a.second < b.second;
                                                 })->first;
        }

        // Odabir najčešćeg znaka za svaku vrstu mutacije
        char base;
        if (votes.substitutionVotes > votes.insertionVotes && votes.substitutionVotes > votes.deletionVotes) {
            base = mostCommonSubstitutionBase;
        } else if (votes.insertionVotes > votes.substitutionVotes && votes.insertionVotes > votes.deletionVotes) {
            base = mostCommonInsertionBase.empty() ? '-' : mostCommonInsertionBase[0];
        } else {
            base = mostCommonDeletionBase.empty() ? '-' : mostCommonDeletionBase[0];
        }

        // Ispis rezultata
        if (base != ' ') {
            if (votes.substitutionVotes >= votes.insertionVotes && votes.substitutionVotes >= votes.deletionVotes) {
                outfile << "substitution,X," << pos << "," << base << "\n";
            } else if (votes.insertionVotes >= votes.substitutionVotes && votes.insertionVotes >= votes.deletionVotes) {
                outfile << "insertion,I," << pos << "," << base << "\n";
            } else {
                outfile << "deletion,D," << pos << "," << base << "\n";
            }
        }
    }

    outfile.close();
}

int main() {
    string samFile = "data/minimap_output.sam";
    string referenceFile = "data/lambda.fasta";
    string outputCsv = "data/mutations.csv";

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

