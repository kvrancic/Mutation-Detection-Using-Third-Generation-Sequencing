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

        if (flag & 4) continue; // preskoči redove gdje je FLAG 4

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
    vector<char> insertionBases;
    int deletionVotes = 0;
    int noneVotes = 0;
};

void addMutationProposal(map<int, MutationProposal>& mutationProposals, int startPos, int readPos, const string& readSeq, const string& reference, char mutationType, int count, ofstream& logFile) {
    logFile << "Adding mutation proposal at position " << startPos << " for type " << mutationType << " with count " << count << endl;
    for (int i = 0; i < count; ++i) {
        int pos = startPos + i;
        if (mutationType == 'M') {
            if (pos < reference.size() && readPos + i < readSeq.size() && reference[pos] == readSeq[readPos + i]) {
                mutationProposals[pos].noneVotes++;
            } else {
                mutationProposals[pos].substitutionVotes++;
                mutationProposals[pos].substitutionBases.push_back(readSeq[readPos + i]);
            }
        } else if (mutationType == 'I') {
            mutationProposals[pos].insertionVotes++;
            if (readPos + i < readSeq.size()) {
                mutationProposals[pos].insertionBases.push_back(readSeq[readPos + i]);
            }
        } else if (mutationType == 'D') {
            mutationProposals[pos].deletionVotes++;
        }
        logFile << "Position: " << pos << ", Mutation Proposal: {substitutionVotes: " << mutationProposals[pos].substitutionVotes
                << ", insertionVotes: " << mutationProposals[pos].insertionVotes << ", deletionVotes: " << mutationProposals[pos].deletionVotes
                << ", noneVotes: " << mutationProposals[pos].noneVotes << "}" << endl;
    }
}

// Funkcija za detekciju mutacija
void detectMutations(const vector<SamEntry>& samEntries, const string& reference, const string& outputCsv, const string& logFilename) {
    ofstream outfile(outputCsv);
    if (!outfile) {
        cerr << "Error opening file: " << outputCsv << endl;
        return;
    }

    ofstream logFile(logFilename);
    if (!logFile) {
        cerr << "Error opening file: " << logFilename << endl;
        return;
    }

    outfile << "type,X,pos,base\n";
    logFile << "type,X,pos,base\n";

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

            // add output to log file so we can see what is happening
            logFile << "--------------------------------------------------------" << endl;
            logFile << "Operation: " << op << ", Count: " << count << endl;
            

            addMutationProposal(mutationProposals, refPos, readPos, readSeq, reference, op, count, logFile);

            if (op == 'M' || op == 'D') {
                refPos += count;
            }
            if (op == 'M' || op == 'I') {
                readPos += count;
            } else if (op == 'S') {
                readPos += count;
            }
        }
    }

    // Većinsko glasanje za svaku poziciju
    for (auto& proposal : mutationProposals) {
        int pos = proposal.first;
        MutationProposal& votes = proposal.second;

        // Log current mutation proposal
        logFile << "--------------------------------------------------------" << endl;
        logFile << "Position: " << pos << ", Mutation Proposal: {substitutionVotes: " << votes.substitutionVotes
                << ", insertionVotes: " << votes.insertionVotes << ", deletionVotes: " << votes.deletionVotes
                << ", noneVotes: " << votes.noneVotes << "}" << endl;

        string type;
        char base = ' ';
        // if none votes is max 
        if (votes.noneVotes >= votes.substitutionVotes && votes.noneVotes >= votes.insertionVotes && votes.noneVotes >= votes.deletionVotes) {
            
            logFile << "None votes won the voting, no mutation detected." << endl;
            
            continue;

        } else if (votes.substitutionVotes >= votes.insertionVotes && votes.substitutionVotes >= votes.deletionVotes && votes.substitutionVotes >= votes.noneVotes) {
            logFile << "Substitution won the voting, possible bases: ";
            for (char base : votes.substitutionBases) {
                logFile << base << " ";
            }
            logFile << endl;

            if (!votes.substitutionBases.empty()) {
                map<char, int> substitutionBaseCounts;
                for (char base : votes.substitutionBases) {
                    substitutionBaseCounts[base]++;
                }
                base = max_element(
                    substitutionBaseCounts.begin(),
                    substitutionBaseCounts.end(),
                    [](const pair<char, int>& a, const pair<char, int>& b) {
                        return a.second < b.second;
                    }
                )->first;
            }
            // inform which base has "won" 
            logFile << "Base " << base << " has won the voting." << endl;

            type = "substitution";
            
            

        } else if (votes.insertionVotes >= votes.substitutionVotes && votes.insertionVotes >= votes.deletionVotes && votes.insertionVotes >= votes.noneVotes) {

            logFile << "Insertion won the voting, possible bases: ";
            for (char base : votes.insertionBases) {
                logFile << base << " ";
            }
            logFile << endl;

            if (!votes.insertionBases.empty()) {
                map<char, int> insertionBaseCounts;
                for (char base : votes.insertionBases) {
                    insertionBaseCounts[base]++;
                }
                base = max_element(
                    insertionBaseCounts.begin(),
                    insertionBaseCounts.end(),
                    [](const pair<char, int>& a, const pair<char, int>& b) {
                        return a.second < b.second;
                    }
                )->first;
            }
            type = "insertion";

            logFile << "Base " << base << " has won the voting." << endl;

        } else if (votes.deletionVotes >= votes.substitutionVotes && votes.deletionVotes >= votes.insertionVotes && votes.deletionVotes >= votes.noneVotes) {
            type = "deletion";
            base = '-';
        } else {
            continue;
        }

        logFile << type << "," << (type == "insertion" ? "I" : (type == "deletion" ? "D" : "X")) << "," << pos << "," << base << "\n";

        // Ispis rezultata
        if (base != ' ') {
            outfile << type << "," << (type == "insertion" ? "I" : (type == "deletion" ? "D" : "X")) << "," << pos << "," << base << "\n";
            
        }
    }

    outfile.close();
    logFile.close();
}

int main() {
    string samFile = "minimap_output.sam";
    string referenceFile = "lambda.fasta";
    string outputCsv = "mutations.csv";
    string logFilename = "mutation_detection_log.txt";

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
    detectMutations(samEntries, reference, outputCsv, logFilename);

    cout << "Mutation detection completed, output saved to " << outputCsv << endl;
    cout << "Log saved to " << logFilename << endl;
    return 0;
}
