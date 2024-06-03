#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>

using std::cerr;
using std::endl;
using std::ifstream;
using std::istringstream;
using std::map;
using std::max_element;
using std::ofstream;
using std::pair;
using std::string;
using std::vector;

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

void addMutationProposal(map<int, MutationProposal>& mutationProposals, int pos, char base, char mutationType, const string& reference, ofstream& logFile) {
    logFile << "Adding mutation proposal at position " << pos << " for type " << mutationType << endl;
    char ref = reference[pos];
    if (mutationType == 'M') {
        if (ref == base) {
            mutationProposals[pos].noneVotes++;
        } else {
            mutationProposals[pos].substitutionVotes++;
            mutationProposals[pos].substitutionBases.push_back(base);
        }
    } else if (mutationType == 'I') {
        mutationProposals[pos].insertionVotes++;
        mutationProposals[pos].insertionBases.push_back(base);
    } else if (mutationType == 'D') {
        mutationProposals[pos].deletionVotes++;
    }
    logFile << "Position: " << pos << ", Mutation Proposal: {substitutionVotes: " << mutationProposals[pos].substitutionVotes
            << ", insertionVotes: " << mutationProposals[pos].insertionVotes << ", deletionVotes: " << mutationProposals[pos].deletionVotes
            << ", noneVotes: " << mutationProposals[pos].noneVotes << "}" << endl;
}

//fja koja gleda koji su mi samentry ukljuceni za tu poziciju, dugo se vrti moozda ako sljedecei nije ukljucena  nije prazna break?
map<pair<int, int>, SamEntry, ComparePositions> getAffectedEntries(const map<pair<int, int>, SamEntry, ComparePositions>& sortedSamEntryMap, int currentPosition) {
    map<pair<int, int>, SamEntry, ComparePositions> affectedEntries;
    for (const auto& entry : sortedSamEntryMap) {
        if (entry.first.first <= currentPosition && currentPosition < entry.first.second) {
            affectedEntries.insert(entry); // Ubacuje par {key, value} u mapu
        }
    }
    return affectedEntries;
}

void adjustSequences(map<pair<int, int>, SamEntry, ComparePositions>& affectedEntries, int currentPosition, char majorityOp, char majorityBase, const string& reference) {
    for (auto& entry : affectedEntries) {
        int refPos = entry.first.first; 
        string& readSeq = entry.second.seq;
        char refBase = reference[currentPosition];

        //int localRefPos = refPos;
        //int readIndex = currentPosition - refPos;
        auto cigarOperations = parseCigar(entry.second.cigar);
        for (const auto& op_count : cigarOperations) {
            int readIndex = currentPosition - entry.first.first;

            char cigarOp = op_count.first;
            int count = op_count.second;
            refPos += count;
            if(refPos >= currentPosition){
                if(cigarOp != 'S'){
                    if (majorityOp == 'I') {
                                if (cigarOp == 'M') {
                                    readSeq.insert(readIndex, 1, majorityBase);
                                    return;
                                } else if (cigarOp == 'D') {
                                    readSeq.insert(readIndex, 1, refBase); //za ovo nisam sigurna
                                    readSeq.insert(readIndex, 1, majorityBase);
                                    return;
                                } else{
                                    return;
                                }
                    } else if (majorityOp == 'D') {
                                if (cigarOp == 'M') {
                                    readSeq.erase(readIndex, 1);
                                    return;
                                } else if (cigarOp == 'I') {
                                    if (readIndex < readSeq.size() - 1) {
                                        readSeq.erase(readIndex, 2);
                                        return;
                                    }
                                } else {
                                    return;
                                }
                    } else if (majorityOp == 'M') {
                                if(majorityBase == ' '){ //tad imam match
                                    majorityBase = refBase;
                                }
                                if (cigarOp == 'I') {
                                    if (readIndex < readSeq.size()) {
                                        readSeq.erase(readIndex, 1);
                                        return;
                                    }
                                } else if (cigarOp == 'D') {
                                    readSeq.insert(readIndex, 1, majorityBase);
                                    return;
                                } else{
                                    return;
                                }
                            }
                } else{
                    return;
                }

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
        //filtriram one zahvacene sam entry
        auto affectedEntries = getAffectedEntries(sortedSamEntryMap, currentPosition);

        if (affectedEntries.empty()) { //ili size manji od 2 pa ako samo jedan onda je on pobjedio za tu poziciju?
            continue;
        }

        for (const auto& entry : affectedEntries) {
            //int refPos = entry.second.pos; // - 1;
            string readSeq = entry.second.seq;

            if (entry.second.flag & 16) {
                std::reverse(readSeq.begin(), readSeq.end());
            }

            auto cigarOperations = parseCigar(entry.second.cigar);
            //int localRefPos = refPos;
            int localReadPos = entry.first.first;
            char op = ' ';
            char base = ' ';

            for (const auto& op_count : cigarOperations) {
                char cigarOp = op_count.first;
                int count = op_count.second;
                localReadPos += count; 
                if(localReadPos >= currentPosition){
                    // presli smo poziciju znaci da je to ta cigar operacija
                    if(cigarOp == 'S'){
                        break;
                    } else {
                    int currentSeq = currentPosition - entry.first.first;
                    base = readSeq[currentSeq];
                    addMutationProposal(mutationProposals, currentPosition, base, cigarOp, reference, logFile);
                    break;
        

                    }
            
        }
        }
        }

        if (mutationProposals[currentPosition].noneVotes > 0 || mutationProposals[currentPosition].substitutionVotes > 0
        || mutationProposals[currentPosition].deletionVotes > 0 || mutationProposals[currentPosition].insertionVotes > 0) {
            const MutationProposal& votes = mutationProposals[currentPosition];
            bool noneAction = false;
            string type;
            char base = ' ';
            if (votes.noneVotes >= votes.substitutionVotes && votes.noneVotes >= votes.insertionVotes && votes.noneVotes >= votes.deletionVotes) {
                //baza je ona na ref
                //base = reference[currentPosition];
                adjustSequences(affectedEntries, currentPosition, 'M', base, reference);
                noneAction = true;
                continue; // nema mutacija
            } else if (votes.substitutionVotes >= votes.insertionVotes && votes.substitutionVotes >= votes.deletionVotes) {
                if (!votes.substitutionBases.empty()) {
                    map<char, int> substitutionBaseCounts;
                    for (char b : votes.substitutionBases) {
                        substitutionBaseCounts[b]++;
                    }
                    base = std::max_element(
                        substitutionBaseCounts.begin(),
                        substitutionBaseCounts.end(),
                        [](const pair<char, int>& a, const pair<char, int>& b) {
                            return a.second < b.second;
                        }
                    )->first;
                }
                type = "substitution";
                adjustSequences(affectedEntries, currentPosition, 'M', base, reference);
            } else if (votes.insertionVotes >= votes.substitutionVotes && votes.insertionVotes >= votes.deletionVotes) {
                if (!votes.insertionBases.empty()) {
                    map<char, int> insertionBaseCounts;
                    for (char b : votes.insertionBases) {
                        insertionBaseCounts[b]++;
                    }
                    base = std::max_element(
                        insertionBaseCounts.begin(),
                        insertionBaseCounts.end(),
                        [](const pair<char, int>& a, const pair<char, int>& b) {
                            return a.second < b.second;
                        }
                    )->first;
                }
                type = "insertion";
                adjustSequences(affectedEntries, currentPosition, 'I', base, reference);
            } else if (votes.deletionVotes >= votes.substitutionVotes && votes.deletionVotes >= votes.insertionVotes) {
                type = "deletion";
                base = '-';
                adjustSequences(affectedEntries, currentPosition, 'D', base, reference);
            } else {
                continue;
            }

            if (base != ' ' && noneAction == false) {
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

   //punim da dobim end position
    map<pair<int, int>, SamEntry> samEntryMap;
    for (const auto& entry : samEntries) {
        int startPos = entry.pos;
        int endPos = startPos;

        auto cigarOperations = parseCigar(entry.cigar);
        for (const auto& op_count : cigarOperations) {
            char op = op_count.first;
            int count = op_count.second;
            endPos += count; 
            
        }

        samEntryMap[{startPos, endPos}] = entry;
    }

    map<pair<int, int>, SamEntry, ComparePositions> sortedSamEntryMap(samEntryMap.begin(), samEntryMap.end());

    ifstream refFile(referenceFile);
    string reference;
    string line;
    while (std::getline(refFile, line)) {
        if (line[0] != '>') {
            reference += line;
        }
    }

    detectMutations(sortedSamEntryMap, reference, outputCsv, logFilename);

    std::cout << "Mutation detection completed, output saved to " << outputCsv << endl;
    std::cout << "Log saved to " << logFilename << endl;
    return 0;
}
