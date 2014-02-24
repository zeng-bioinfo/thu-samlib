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
#include <map>
using namespace std;


namespace HomopolymerSpace{
    const int HomopolymerSizeMax=100;
    const int HomopolymerPosMax=20;
    using namespace NucleotideSpace;
    const int ALPHABITSIZE=3;
    const int BETABITSIZE=2;
    const int LENGTHBITSIZE=7;
    const int QUALITYBITSIZE=6;
    const int STRANDBITSIZE=1;

    struct HomopolymerSequence{ // Homopolymer sequence object
        vector<string> alpha;   // nucleotide of the homopolymer segment
        vector<int> ell;        // length of the homopolymer segment
        vector<int> t0;         // starting position
        vector<int> t1;         // ending position
        int len;                // length of the homopolymer sequence

        // Default constructor
        HomopolymerSequence(){
            alpha=vector<string>();
            ell=vector<int>();
            t0=vector<int>();
            t1=vector<int>();
            len=0;
        }

        // Constructor from a given nucleotide sequence
        HomopolymerSequence(const string& nucl_seq){
            alpha=vector<string>();
            ell=vector<int>();
            t0=vector<int>();
            t1=vector<int>();
            len=0;
            homopolymerize(nucl_seq);
        }

        // Default destructor
        virtual ~HomopolymerSequence(){
            alpha.clear();
            ell.clear();
            t0.clear();
            t1.clear();
        }

        // Transformer to convert a nucleotide sequence to a homopolymer sequence
        void homopolymerize(const string& nucl_seq){
            string curr_alpha="";   // alpha of current homopolymer
            int curr_ell=0;         // length of current homopolymer
            int curr_t0=0, curr_t1=-1;

            for (int i=0; i<(int)nucl_seq.length(); i++){
                // update current homopolymer information
                curr_alpha=nucl_seq.substr(i,1);
                curr_ell+=1;
                curr_t1+=1;
                // check whether move on to next position
                if (i+1>=(int)nucl_seq.length() ||
                        (curr_alpha!=nucl_seq.substr(i+1,1) && i+1<(int)nucl_seq.length())){
                    alpha.push_back(curr_alpha);
                    ell.push_back(curr_ell);
                    t0.push_back(curr_t0);
                    t1.push_back(curr_t1);
                    len+=1;

                    // fresh the memory
                    curr_alpha="";
                    curr_ell=0;
                    curr_t0=i+1;
                    curr_t1=i;
                }
            }
        }
    };
}


namespace SemiHomopolymerAlignmentSpace{
    enum AlignmentState{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,A18,A19,A20,    // g15 is "greater than 15"
                        C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,    // g15 is "greater than 15"
                        G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,G16,G17,G18,G19,G20,    // g15 is "greater than 15"
                        T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,T16,T17,T18,T19,T20,    // g15 is "greater than 15"
                        I,SH,ST,B,E};
    const int ALIGNMENTSTATSIZE=4*20+1+2+2;
    const int HOMOPOLYMERSIZEMAX=100;

    extern std::map<string, int> state2idx;
    extern std::map<int, string> idx2state;
    extern std::map<int, int> state_length;
    extern std::map<int, char> state_alpha;

    /**
     * @brief The Query struct
     */
    struct Query{                   // Query is the nucleotide sequence
        vector<string> seq;
        vector<string> qual;
        int len;
        Query():len(0){}
        Query& operator=(const Query& q){
            reset();
            seq=q.seq;
            qual=q.qual;
            len=q.len;
            return *this;
        }
        void reset(){
            vector<string>().swap(seq);
            vector<string>().swap(qual);
            len=0;
        }
        void push_back(string _s,string _q){
            seq.push_back(_s);
            qual.push_back(_q);
            len+=1;
        }
        void insert_head(string _s, string _q){
            seq.insert(seq.begin(),_s);
            qual.insert(qual.begin(),_q);
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
            reset();
            alpha=t.alpha;
            ell=t.ell;
            len=t.len;
            return *this;
        }
        void reset(){
            vector<string>().swap(alpha);
            vector<int>().swap(ell);
            len=0;
        }
        void push_back(string _a,int _l){
            alpha.push_back(_a);
            ell.push_back(_l);
            len+=1;
        }
        void insert_head(string _a, int _l){
            alpha.insert(alpha.begin(),_a);
            ell.insert(ell.begin(), _l);
            len+=1;
        }
    };
    /**
     * @brief The Status struct
     */
    struct Status{                  // Status is the status sequence
        vector<AlignmentState> status;
        int len;
        Status():len(0){}
        Status& operator=(const Status& s){
            reset();
            status=s.status;
            len=s.len;
            return *this;
        }
        void reset(){
            vector<AlignmentState>().swap(status);
            len=0;
        }
        void push_back(AlignmentState _s){
            status.push_back(_s);
            len+=1;
        }
        void insert_head(AlignmentState _s){
            status.insert(status.begin(),_s);
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
    // TODO: higher order markov model
    // TODO: compute new quality statistics, depending on
    //       match or mismatch
    //       frequency( quality | pi )-->frequency( quality | beta,pi )

    // LOG:
    // @2014.1.27 change the counting table from <alpha,ell>-based to state(pi)-based
    // @2014.2.14 add a new table to count the occurrence of events based on the position on homopolymer and quality score

    public:
        vector<int> base_call_table;    // NEW: frequency( beta | pi ) //OLD:frequency( beta | alpha,l )
        vector<int> qual_call_table;    // NEW: frequency( quality | beta,pi ) //OLD:frequency( quality | alpha,l )
        vector<int> len_call_table;     // NEW: frequency( length | pi ) //OLD:frequency( k | alpha,l )
        vector<int> pos_qual_base_call_table;   // NEW NEW: frequency( alpha->beta | pos, qual, pi )

    public:
        vector<int> state_trans_table;  // frequency( pi | ppi )

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

        // Old routines to count the occurrences of patterns among the alignment
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

        // New routines to count the occurrences of the patterns among the alignment
    public:
        /**
         * @brief bc_incr1
         * @param strand
         * @param pi
         * @param beta
         */
        void bc_incr1(int strand, int pi, int pos, int qual, char beta);
        /**
         * @brief lc_incr1
         * @param strand
         * @param pi
         * @param kappa
         */
        void lc_incr1(int strand, int pi, int kappa);
        /**
         * @brief qc_incr1
         * @param strand
         * @param pi
         * @param qual
         */
        void qc_incr1(int strand, int pi, int pos, int qual);
        /**
         * @brief pqbc_incr1
         * @param strand
         * @param pi
         * @param pos
         * @param qual
         * @param beta
         */
        void pqbc_incr1(int strand, int pi, int pos, int qual, int beta);
        /**
         * @brief bc_elem1
         * @param strand
         * @param pi
         * @param beta
         * @return
         */
        int bc_elem1(int strand, int pi, int pos, int qual, char beta);
        /**
         * @brief lc_elem1
         * @param strand
         * @param pi
         * @param kappa
         * @return
         */
        int lc_elem1(int strand, int pi, int kappa);
        /**
         * @brief qc_elem1
         * @param strand
         * @param pi
         * @param qual
         * @return
         */
        int qc_elem1(int strand, int pi, int pos, int qual);
        /**
         * @brief pqbc_elem1
         * @param strand
         * @param pi
         * @param pos
         * @param qual
         * @param beta
         * @return
         */
        int pqbc_elem1(int strand, int pi, int pos, int qual, int beta);

    public:
        /**
         * @brief st_incr1
         * @param pi
         * @param ppi
         */
        void st_incr1(int strand, int pi, int ppi);
        /**
         * @brief st_elem1
         * @param pi
         * @param ppi
         * @return
         */
        int st_elem1(int strand, int pi, int ppi);

    public:
        /**
         * @brief print
         * @param c
         */
        void print(int c);
        /**
         * @brief print
         * @param filename
         */
        void print(ofstream& ofs, int c);

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
        /**
         * @brief print
         * @param filename
         */
        void print(string filename);

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
        /**
         * @brief open
         * @param filename
         * @param total
         */
        void open(string filename, int total);

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
        /**
         * @brief print
         * @param filename
         */
        void print(string filename);

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
