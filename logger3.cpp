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

vector<pair<char, int>> parseCigar(const string& cigar) {
    vector<pair<char, int>> operations;
    istringstream iss(cigar);
    char op;
    int count;
    while (iss >> count >> op) {
        operations.push_back(make_pair(op, count));
    }
    return operations;
}

struct ComparePositions {
    bool operator()(const pair<int, int>& a, const pair<int, int>& b) const {
        if (a.first != b.first) {
            return a.first < b.first; // sort po pocetnoj poziciji
        } else {
            return a.second < b.second; // jednake, sortiraj po zavrsnoj poziciji
        }
    }
};

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

vector<SamEntry> getAffectedEntries(const map<pair<int, int>, SamEntry, ComparePositions>& sortedSamEntryMap, int currentPosition) {
    vector<SamEntry> affectedEntries;
    for (auto it = sortedSamEntryMap.lower_bound({currentPosition, 0}); it != sortedSamEntryMap.end() && currentPosition < it->first.second; ++it) {
        affectedEntries.push_back(it->second);
    }
    return affectedEntries;
}

void adjustSequences(vector<SamEntry>& affectedEntries, int currentPosition, char majorityOp, char majorityBase, const string& reference) {
    for (auto& entry : affectedEntries) {
        int refPos = entry.pos - 1;
        if (currentPosition < refPos || currentPosition >= refPos + entry.seq.size()) continue;

        string& readSeq = entry.seq;
        int readPos = currentPosition - refPos;

        if (majorityOp == 'I') {
            if (readPos < readSeq.size()) {
                readSeq.insert(readPos, 1, majorityBase);
            }
        } else if (majorityOp == 'D') {
            if (readPos < readSeq.size()) {
                readSeq.erase(readPos, 1);
            }
        } else if (majorityOp == 'M') {
            if (readPos < readSeq.size()) {
                readSeq[readPos] = majorityBase;
            }
        }
    }
}

void detectMutations(map<pair<int, int>, SamEntry, ComparePositions>& sortedSamEntryMap, const string& reference, const string& outputCsv, const string& logFilename) {
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

    for (int currentPosition = 0; currentPosition < reference.size(); ++currentPosition) {
        map<int, MutationProposal> mutationProposals;
        auto affectedEntries = getAffectedEntries(sortedSamEntryMap, currentPosition);

        for (const auto& entry : affectedEntries) {
            int refPos = entry.pos - 1;
            int readPos = currentPosition - refPos;
            string readSeq = entry.seq;

            if (entry.flag & 16) {
                reverse(readSeq.begin(), readSeq.end());
            }

            auto cigarOperations = parseCigar(entry.cigar);
            int localRefPos = refPos;
            int localReadPos = 0;
            for (const auto& op_count : cigarOperations) {
                char op = op_count.first;
                int count = op_count.second;

                if (currentPosition >= localRefPos && currentPosition < localRefPos + count) {
                    addMutationProposal(mutationProposals, localRefPos, localReadPos, readSeq, reference, op, currentPosition - localRefPos + 1, logFile);
                    break;
                }

                if (op == 'M' || op == 'D') {
                    localRefPos += count;
                }
                if (op == 'M' || op == 'I') {
                    localReadPos += count;
                } else if (op == 'S') {
                    localReadPos += count;
                }
            }
        }

        if (mutationProposals.find(currentPosition) != mutationProposals.end()) {
            const MutationProposal& votes = mutationProposals[currentPosition];

            string type;
            char base = ' ';
            if (votes.noneVotes >= votes.substitutionVotes && votes.noneVotes >= votes.insertionVotes && votes.noneVotes >= votes.deletionVotes) {
                continue; // nema mutacija
            } else if (votes.substitutionVotes >= votes.insertionVotes && votes.substitutionVotes >= votes.deletionVotes && votes.substitutionVotes >= votes.noneVotes) {
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
                type = "substitution";
            } else if (votes.insertionVotes >= votes.substitutionVotes && votes.insertionVotes >= votes.deletionVotes && votes.insertionVotes >= votes.noneVotes) {
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
                adjustSequences(affectedEntries, currentPosition, 'I', base, reference);
            } else if (votes.deletionVotes >= votes.substitutionVotes && votes.deletionVotes >= votes.insertionVotes && votes.deletionVotes >= votes.noneVotes) {
                type = "deletion";
                base = '-';
                adjustSequences(affectedEntries, currentPosition, 'D', base, reference);
            } else {
                continue;
            }

            if (base != ' ') {
                outfile << type << "," << (type == "insertion" ? "I" : (type == "deletion" ? "D" : "X")) << "," << currentPosition << "," << base << "\n";
            }
        }
    }

    outfile.close();
    logFile.close();
}

int main() {
    string samFile = "data/minimap_output.sam";
    string referenceFile = "data/lambda.fasta";
    string outputCsv = "data/mutations.csv";
    string logFilename = "data/mutation_detection_log.txt";

    vector<SamEntry> samEntries = parseSamFile(samFile);

    map<pair<int, int>, SamEntry> samEntryMap;
    for (const auto& entry : samEntries) {
        int startPos = entry.pos;
        int endPos = startPos;

        auto cigarOperations = parseCigar(entry.cigar);
        for (const auto& op_count : cigarOperations) {
            char op = op_count.first;
            int count = op_count.second;
            if (op == 'M' || op == 'D' || op == 'N' || op == '=' || op == 'X') {
                endPos += count;
            }
        }

        samEntryMap[{startPos, endPos}] = entry;
    }

    map<pair<int, int>, SamEntry, ComparePositions> sortedSamEntryMap(samEntryMap.begin(), samEntryMap.end());

    ifstream refFile(referenceFile);
    string reference;
    string line;
    while (getline(refFile, line)) {
        if (line[0] != '>') {
            reference += line;
        }
    }

    detectMutations(sortedSamEntryMap, reference, outputCsv, logFilename);

    cout << "Mutation detection completed, output saved to " << outputCsv << endl;
    cout << "Log saved to " << logFilename << endl;
    return 0;
}
