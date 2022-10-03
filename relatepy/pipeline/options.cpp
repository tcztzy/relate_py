#include "cxxopts.hpp"

// I found it is difficult for me to directly call following codes, so I write this.

cxxopts::Options get_options(int argc, char* argv[]) {
    cxxopts::Options options("Relate");
    options.add_options()
        ("help", "Print help.")
        ("mode", "Choose which part of the algorithm to run.", cxxopts::value<std::string>()) 
        ("haps", "Filename of haps file (Output file format of Shapeit).", cxxopts::value<std::string>())
        ("sample", "Filename of sample file (Output file format of Shapeit).", cxxopts::value<std::string>())
        ("map", "Genetic map.", cxxopts::value<std::string>())
        ("m,mutation_rate", "Mutation rate.", cxxopts::value<float>())
        ("N,effectiveN", "Effective population size.", cxxopts::value<float>())
        ("o,output", "Filename of output without file extension.", cxxopts::value<std::string>())
        ("dist", "Optional but recommended. Distance in BP between SNPs. Can be generated using RelateFileFormats. If unspecified, distances in haps are used.", cxxopts::value<std::string>())
        ("annot", "Optional. Filename of file containing additional annotation of snps. Can be generated using RelateFileFormats.", cxxopts::value<std::string>()) 
        ("memory", "Optional. Approximate memory allowance in GB for storing distance matrices. Default is 5GB.", cxxopts::value<float>())
        ("sample_ages", "Optional. Filename of file containing sample ages (one per line).", cxxopts::value<std::string>()) 
        ("chunk_index", "Optional. Index of chunk. (Use when running parts of the algorithm on an individual chunk.)", cxxopts::value<int>())
        ("first_section", "Optional. Index of first section to infer. (Use when running parts of algorithm on an individual chunk.)", cxxopts::value<int>())
        ("last_section", "Optional. Index of last section to infer. (Use when running parts of algorithm on an individual chunk.)", cxxopts::value<int>())
        ("coal", "Optional. Filename of file containing coalescent rates. If specified, it will overwrite --effectiveN.", cxxopts::value<std::string>()) 
            ("fb", "Optional. Force build a new tree every x bases.", cxxopts::value<float>()) 
        //("anc_allele_unknown", "Specify if ancestral allele is unknown.") 
            ("transversion", "Only use transversion for bl estimation.")
        ("i,input", "Filename of input.", cxxopts::value<std::string>())
            ("painting", "Optional. Copying and transition parameters in chromosome painting algorithm. Format: theta,rho. Default: 0.025,1.", cxxopts::value<std::string>())
        ("seed", "Optional. Seed for MCMC in branch lengths estimation.", cxxopts::value<int>());
    options.parse(argc, argv);
    return options;
}
