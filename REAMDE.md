# Mutation Detection Using Third-Generation Sequencing

## Project Overview
This is a university project for a bioinformatics course focused on detecting mutations in DNA sequences using third-generation sequencing data. The goal is to map sequencing reads to a reference genome, align them, and identify various types of mutations, including substitutions, insertions, and deletions. This project is implemented in C++.

## Installation

### Prerequisites
- C++ Compiler (GCC recommended)
- Minimap2
- FreeBayes

### Installing Minimap2
Clone the Minimap2 repository and compile it:

    git clone https://github.com/lh3/minimap2
    cd minimap2 && make

### Installing FreeBayes
Clone the FreeBayes repository and compile it:

    git clone --recursive https://github.com/freebayes/freebayes.git
    cd freebayes
    make

## Usage
Instructions on how to run the project with an example command:

    # Example command
    ./mutation_detector --reference genome.fa --reads reads.fa --output mutations.csv

## Input and Output Formats

### Input
- **Reference Genome (FASTA format)**: Reference genome against which reads will be mapped.
- **Mutated Genome Reads (FASTA format)**: Sequencing reads from the mutated genome.

### Output
- **CSV File**: A CSV file listing the detected mutations with details like position and type of mutation.

## Example Output
Example of what the output might look like:

    Mutation, Position, Base
    Substitution, 100, A>T
    Insertion, 150, -C
    Deletion, 200, -

## Contributing
This project is part of a bioinformatics course and is primarily for educational purposes. Contributions and suggestions are welcome to improve the project or address issues.

## License
Specify the license under which your project is made available.

## Authors
- Tea Ćetojević-Tisaj 
- Karlo Vrančić

## Acknowledgments
- This project is developed as a part of the coursework for the Bioinformatics course at FER.
- Thanks to the provided resources and tools like Minimap2 and FreeBayes that are integral to the project's implementation.
