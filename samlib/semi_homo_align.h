#ifndef SEMIHOMOPOLYMERALIGNMENT_H
#define SEMIHOMOPOLYMERALIGNMENT_H

#include "nucl_align.h"
#include <string>
#include <vector>

namespace SemiHomopolymerSpace{
    using namespace std;
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
/**
 * @brief The SemiHomopolymerAlignment class
 */
using namespace std;
using namespace SemiHomopolymerSpace;
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
