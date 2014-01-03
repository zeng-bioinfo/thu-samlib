#ifndef NUCL_ALIGN_H
#define NUCL_ALIGN_H

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;

/**
 *
 */
namespace NucleotideSpace {
    /**
    * @brief The Nucleotide enum
    */
    enum Nucleotide{A,C,G,T,S};
    /**
     * @brief The Metrics enum
     */
    enum Metrics{NucleotideSize=5,QualityScoreMax=60};
}

/**
 * @brief The NucleotideAlignmentStat class
 *        It is to store the statistics of an alignment
 */
class NucleotideAlignmentStat{
    // 0-order statistics
    public:
        /**
         * @brief base_call_table
         *
         *   OUTPUT    | A | C | G | T | - |
         *----------------------------------
         *   I     | A |0,0|0,1|0,2|0,3|0,4|
         *   N     | C |1,0|1,1|1,2|1,3|1,4|
         *   P     | G |2,0|2,1|2,2|2,3|2,4|
         *   U     | T |3,0|3,1|3,2|3,3|3,4|
         *   T     | - |4,0|4,1|4,2|4,3|4,4|
         *
         */
        int base_call_table[NucleotideSpace::NucleotideSize][NucleotideSpace::NucleotideSize];

        /**
         * @brief qual_call_table
         */
        int qual_call_table[NucleotideSpace::NucleotideSize][NucleotideSpace::QualityScoreMax];

    public:
        NucleotideAlignmentStat& operator=(const NucleotideAlignmentStat& _align_stat);

    public:
        NucleotideAlignmentStat();  // default constructor
        NucleotideAlignmentStat(const NucleotideAlignmentStat& _align_stat);    // copy constructor
        virtual ~NucleotideAlignmentStat(); // default deconstructor
};

/**
 * @brief The NucleotideAlignment class
 *        It is to represent the alignment of query to target in nucleotide space.
 */
class NucleotideAlignment
{
    // variable member
    public:
        string align_name;  // alignment id
        string query_seq;   // nucleotide sequence
        string query_qual;  // base quality score
        string target_seq;  // nucleotide sequence
        string align_status;  // status in every aligned position

    // function member
    public:
        // overload operator=
        NucleotideAlignment& operator=(const NucleotideAlignment& _align);

    public:
        /**
         * @brief indel_shift_right
         *        Shift the space symbols to right end as possible
         */
        void indel_shift_right();

    public:
        /**
         * @brief statistics
         *        Print out the statistics of the alignment
         */
        void statistics();
        /**
         * @brief statistics
         *        Get the statistics of the alignment, ignore the sequencing cycle
         * @param stat
         */
        void statistics(NucleotideAlignmentStat& stat);
        /**
         * @brief statistics
         *        Get the statistics of the alignment, consider the sequencing cycle
         * @param stat
         * @param cycles
         */
        void statistics(vector<NucleotideAlignmentStat>& stat, int cycles);



    public:
        NucleotideAlignment();  // default constructor
        NucleotideAlignment(const NucleotideAlignment& _align); // copy constructor
        virtual ~NucleotideAlignment(); // default deconstructor
};


/**
 * @brief The NucleotideAlignmentPool class
 */
class NucleotideAlignmentPool{
    // a set of NucleotideAlignments
    public:
        vector<NucleotideAlignment> align_pool;
        NucleotideAlignmentStat align_stat;

    public:
        /**
         * @brief open
         *        To load in the alignments
         * @param filename
         */
        void open(string filename);

    public:
        NucleotideAlignmentPool();  // default constructor
        ~NucleotideAlignmentPool(); // default deconstructor
};

#endif // NUCL_ALIGN_H
