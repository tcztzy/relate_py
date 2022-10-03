cdef extern from "cxxopts.hpp" namespace "cxxopts":
    cdef cppclass Options:
        pass

cdef extern from "options.cpp":
    Options get_options(int argc, char** argv)

cdef extern from "Paint.cpp":
    int Paint(Options& options, int chunk_index)

cdef extern from "BuildTopology.cpp":
    int BuildTopology(Options& options, int chunk_index, int first_section, int last_section)

cdef extern from "FindEquivalentBranches.cpp":
    int FindEquivalentBranches(Options& options, int chunk_index)

cdef extern from "InferBranchLengths.cpp":
    int GetBranchLengths(Options& options, int chunk_index, int first_section, int last_section)

cdef extern from "CombineSections.cpp":
    int CombineSections(Options& options, int chunk_index)
    int CombineSections(Options& options)

cdef extern from "Finalize.cpp":
    int Finalize(Options&)
