# cython: language_level=3
from libc.stdio cimport FILE
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "collapsed_matrix.hpp":
    cdef cppclass CollapsedMatrix[T]:
        pass

cdef extern from "data.hpp":
    cdef cppclass Data:
        string name
        int N  # number of sequences
        int L  # number of SNPs
        int Ne  # effective population size
        double mu  # mutation rate
        double theta  # mutation probability for painting. set to 0.001
        CollapsedMatrix[char] sequence  # sequence matrix, containing 0 and 1
        vector[double] r  # vector of recombination distances from one SNP to the next
        Data(const char* filename_sequence, const char* filename_pos, const char* filename_dist, const char* filename_rec, const char* filename_rpos, const char* filename_state, int Ne, double mu) except +
        Data(const char* filename_sequence, const char* filename_pos, const char* filename_dist, const char* filename_rec, const char* filename_rpos, const char* filename_state, int Ne) except +
        Data(const char* filename_sequence, const char* filename_pos, const char* filename_dist, const char* filename_rec, const char* filename_rpos, const char* filename_state) except +
        Data(const char* filename_dist, const char* filename_param, int Ne, double mu) except +
        Data(const char* filename_dist, const char* filename_param, int Ne) except +
        Data(const char* filename_dist, const char* filename_param) except +
        Data() except +
        void MakeChunks(const string& filename_haps, const string& filename_sample, const string& filename_map, const string& filename_dist, const string& file_out, bool use_transition, float max_memory)
        void MakeChunks(const string& filename_haps, const string& filename_sample, const string& filename_map, const string& filename_dist, const string& file_out, bool use_transition)
