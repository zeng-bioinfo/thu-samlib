// TODO List
// [1]  Calculate the statistics depending on the alignment status
// [2]  Calculate the statistics depending on the read strand           [Done]
// [3]  Calculate the higher order markov statistics
// ...

#ifndef SEMIHOMOPOLYMERALIGNMENT_H
#define SEMIHOMOPOLYMERALIGNMENT_H

#include "nucl_align.h"
#include <string>
#include <vector>
using namespace std;

namespace HomopolymerSpace{
    const int HomopolymerSizeMax=100;
    using namespace NucleotideSpace;
    const int ALPHABITSIZE=3;
    const int BETABITSIZE=2;
    const int LENGTHBITSIZE=7;
    const int QUALITYBITSIZE=6;
    const int STRANDBITSIZE=1;
}

namespace SemiHomopolymerAlignmentSpace{
    /**
     * @brief The Query struct
     */
    struct Query{                   // Query is the nucleotide sequence
        vector<string> seq;
        vector<string> qual;
        int len;
        Query():len(0){}
        Query& operator=(const Query& q){
            seq=q.seq;
            qual=q.qual;
            len=q.len;
            return *this;
        }
        void reset(){
            seq=vector<string>();
            qual=vector<string>();
            len=0;
        }
        void push_bach(string _s,string _q){
            seq.push_back(_s);
            qual.push_back(_q);
            len+=1;
        }
    };
    /**
     * @brief The Target struct
     */
    struct Target{                  // Target is the homopolymer sequence
        vector<string> alpha;       // alpha is the nucleotide of a homopolymer
        vector<int> ell;            // ell is the length of a homopolymer
        int len;
        Target():len(0){}
        Target& operator=(const Target& t){
            alpha=t.alpha;
            ell=t.ell;
            len=t.len;
            return *this;
        }
        void reset(){
            alpha=vector<string>();
            ell=vector<int>();
            len=0;
        }
        void push_back(string _a,int _l){
            alpha.push_back(_a);
            ell.push_back(_l);
            len+=1;
        }
    };
    /**
     * @brief The Status struct
     */
    struct Status{                  // Status is the status sequence
        vector<string> status;
        int len;
        Status():len(0){}
        Status& operator=(const Status& s){
            status=s.status;
            len=s.len;
            return *this;
        }
        void reset(){
            status=vector<string>();
            len=0;
        }
        void push_back(string _s){
            status.push_back(_s);
            len+=1;
        }
    };
}
using namespace SemiHomopolymerAlignmentSpace;

/**
 * @brief The SemiHomopolymerAlignmentStat class
 */
class SemiHomopolymerAlignmentStat
{
    // TODO: consider the alignment status
    // TODO: higher order markov model

    public:
        vector<int> base_call_table;    // frequency( beta | alpha,l )
        vector<int> qual_call_table;    // frequency( quality | alpha,l )
        vector<int> len_call_table;     // frequency( k | alpha,l )

    public:
        /**
         * @brief operator =
         * @param another
         * @return
         */
        SemiHomopolymerAlignmentStat& operator=(const SemiHomopolymerAlignmentStat& another);
        /**
         * @brief operator +=
         * @param another
         * @return
         */
        SemiHomopolymerAlignmentStat& operator+=(const SemiHomopolymerAlignmentStat& another);

    public:
        /**
         * @brief bc_incr1
         *        increase an element of base-call table by one
         * @param alpha
         * @param ell
         * @param beta
         */
        void bc_incr1(int strand, char alpha, int ell, char beta);
        /**
         * @brief qc_incr1
         *        increase an element of quality-score-call table by one
         * @param alpha
         * @param ell
         * @param qual
         */
        void qc_incr1(int strand, char alpha, int ell, int qual);
        /**
         * @brief lc_incr1
         *        increase an element of length-call table by one
         * @param alpha
         * @param ell
         * @param kappa
         */
        void lc_incr1(int strand, char alpha, int ell, int kappa);

        /**
         * @brief bc_elem1
         * @param alpha
         * @param ell
         * @param beta
         * @return
         */
        int bc_elem1(int strand, char alpha, int ell, char beta);
        /**
         * @brief qc_elem1
         * @param alpha
         * @param ell
         * @param qual
         * @return
         */
        int qc_elem1(int strand, char alpha, int ell, int qual);
        /**
         * @brief lc_elem1
         * @param alpha
         * @param ell
         * @param kappa
         * @return
         */
        int lc_elem1(int strand, char alpha, int ell, int kappa);

    public:
        /**
         * @brief print
         */
        void print();
        /**
         * @brief print
         * @param filename
         */
        void print(string filename);

    // constructor and destructor
    public:
        /**
        * @brief SemiHomopolymerAlignmentStat
        */
        SemiHomopolymerAlignmentStat();
        /**
        * @brief SemiHomopolymerAlignmentStat
        * @param another
        */
        SemiHomopolymerAlignmentStat(const SemiHomopolymerAlignmentStat& another);
        /**
        * @brief ~SemiHomopolymerAlignmentStat
        */
        virtual ~SemiHomopolymerAlignmentStat();
};

/**
 * @brief The SemiHomopolymerAlignment class
 */
class SemiHomopolymerAlignment:public NucleotideAlignment
{
    public:
        Query align_query;
        Target align_target;
        Status align_status;

    // conversion
    public:
        /**
         * @brief homopolymerize
         * @param nucl_align
         */
        void homopolymerize(const NucleotideAlignment& nucl_align);

    // overloading operators
    public:
        /**
         * @brief operator =
         * @param homo_align
         * @return
         */
        SemiHomopolymerAlignment& operator=(const SemiHomopolymerAlignment& homo_align);

        /**
         * @brief operator =
         * @param nucl_align
         * @return
         */
        SemiHomopolymerAlignment& operator=(const NucleotideAlignment& nucl_align);

    public:
        /**
         * @brief print
         */
        void print();

    public:
        /**
         * @brief statistics
         * @param stat
         * @param cycles
         */
        void statistics(vector<SemiHomopolymerAlignmentStat>& stat, int cycles);

    // constructor & destructor
    public:
        /**
         * @brief SemiHomopolymerAlignment
         */
        SemiHomopolymerAlignment();
        /**
         * @brief SemiHomopolymerAlignment
         * @param homo_align
         */
        SemiHomopolymerAlignment(const SemiHomopolymerAlignment& homo_align);
        /**
         * @brief SemiHomopolymerAlignment
         * @param nucl_align
         */
        SemiHomopolymerAlignment(const NucleotideAlignment& nucl_align);
        /**
        * @brief ~SemiHomopolymerAlignment
        */
        virtual ~SemiHomopolymerAlignment();
};

class SemiHomopolymerAlignmentPool{
    public:
        /**
         * @brief align_pool
         */
        vector<SemiHomopolymerAlignment> align_pool;

    public:
        /**
         * @brief open
         * @param filename
         */
        void open(string filename);

    public:
        /**
         * @brief operator =
         * @param another_pool
         * @return
         */
        SemiHomopolymerAlignmentPool& operator=(const SemiHomopolymerAlignmentPool& another);

    public:
        /**
         * @brief print
         */
        void print();

    public:
        /**
         * @brief statistics
         * @param cycles
         */
        void statistics(int cycles);
        /**
         * @brief statistics
         * @param filename
         * @param cycles
         */
        void statistics(string filename, int cycles);
        /**
         * @brief statistics
         * @param stat
         * @param cycles
         */
        void statistics(vector<SemiHomopolymerAlignmentStat>& stat, int cycles);

    public:
        /**
         * @brief SemiHomopolymerAlignment
         */
        SemiHomopolymerAlignmentPool();

        /**
         * @brief SemiHomopolymerAlignmentPool
         * @param another_pool
         */
        SemiHomopolymerAlignmentPool(const SemiHomopolymerAlignmentPool& another);

        /**
         * @brief ~SemiHomopolymerAlignmentPool
         */
        virtual ~SemiHomopolymerAlignmentPool();

};

#endif // SEMIHOMOPOLYMERALIGNMENT_H
