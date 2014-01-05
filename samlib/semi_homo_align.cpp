#include "semi_homo_align.h"
#include <iostream>
using namespace std;

/**
 * @brief SemiHomopolymerAlignment::SemiHomopolymerAlignment
 */
SemiHomopolymerAlignment::SemiHomopolymerAlignment(){

}

/**
 * @brief SemiHomopolymerAlignment::SemiHomopolymerAlignment
 * @param _align
 */
SemiHomopolymerAlignment::SemiHomopolymerAlignment(const SemiHomopolymerAlignment &homo_align)
{
    align_name=homo_align.align_name;
    align_query=homo_align.align_query;
    align_target=homo_align.align_target;
    align_status=homo_align.align_status;
}

/**
 * @brief SemiHomopolymerAlignment::SemiHomopolymerAlignment
 * @param nucl_align
 */
SemiHomopolymerAlignment::SemiHomopolymerAlignment(const NucleotideAlignment &nucl_align){
    this->homopolymerize(nucl_align);
}

/**
 * @brief SemiHomopolymerAlignment::~SemiHomopolymerAlignment
 */
SemiHomopolymerAlignment::~SemiHomopolymerAlignment(){
    ;
}

/**
 * @brief SemiHomopolymerAlignment::operator =
 * @param _align
 * @return
 */
SemiHomopolymerAlignment& SemiHomopolymerAlignment::operator =(const SemiHomopolymerAlignment& homo_align){
    align_name=homo_align.align_name;
    align_query=homo_align.align_query;
    align_target=homo_align.align_target;
    align_status=homo_align.align_status;
    return *this;
}

SemiHomopolymerAlignment& SemiHomopolymerAlignment::operator =(const NucleotideAlignment& nucl_align){
    this->homopolymerize(nucl_align);
    return *this;
}

/**
 * @brief SemiHomopolymerAlignment::homopolymerize
 * @param nucl_align
 */
void SemiHomopolymerAlignment::homopolymerize(const NucleotideAlignment &nucl_align){
    // TODO: consider the alignment status
    // TODO: consider IUPAC code?

    // clear everything
    this->align_name="";
    this->align_query.reset();
    this->align_target.reset();
    this->align_status.reset();


    string qis; // query subsequence at the current position
    string qiq; // quality scores of query subsequence

    string tia; // alpha of target subsequence (homopolymer) at the current position
    int til;    // length of target subsequence (homopolymer) at the current position

    // initialization
    qis=""; qiq="";
    tia=""; til=0;

    // scan over the nucleotide alignment
    this->align_name=nucl_align.align_name;
    for (int i=0; i<(int)nucl_align.align_status.length(); i++){
        // update query subsequence
        if (nucl_align.query_seq[i]!='-'){
            qis+=nucl_align.query_seq[i];
            qiq+=nucl_align.query_qual[i];
        }
        // update target subsequence
        if (nucl_align.target_seq[i]!='-'){
            tia=nucl_align.target_seq[i];
            til+=1;
        }
        // check whether the next position is not feasible
        if (i+1>=(int)nucl_align.align_status.length() ||
                (nucl_align.target_seq[i+1]!='-' && tia!=nucl_align.target_seq.substr(i+1,1) && i+1<(int)nucl_align.align_status.length()) ||
                (nucl_align.target_seq[i+1]=='-' && tia!=nucl_align.query_seq.substr(i+1,1) && i+1<(int)nucl_align.align_status.length())){
            // package up a complete homopolymer position
            this->align_query.push_bach(qis,qiq);
            this->align_target.push_back(tia,til);
            // reset temporaty buffer
            qis=""; qiq="";
            tia=""; til=0;
        }
    }
}

/**
 * @brief SemiHomopolymerAlignment::print
 */
void SemiHomopolymerAlignment::print(){
    // TODO: consider the alignment status
    // ...

    // temporary buffer
    string buffer;

    // print out message
    cout<<this->align_name<<endl;

    // print out query sequence
    for (int i=0; i<align_query.len-1; i++){
        buffer=align_query.seq[i];
        if (align_target.ell[i]>(int)buffer.length()){
            for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                buffer+=" ";
            }
        }
        cout<<buffer<<"  ";
    }
    buffer=align_query.seq[align_query.len-1];
    if (align_target.ell[align_target.len-1]>(int)buffer.length()){
        for (int j=align_target.ell[align_query.len-1]-buffer.length(); j>0; j--){
            buffer+=" ";
        }
    }
    cout<<buffer<<endl;

    // print out query quality score
    for (int i=0; i<align_query.len-1; i++){
        buffer=align_query.qual[i];
        if (align_target.ell[i]>(int)buffer.length()){
            for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                buffer+=" ";
            }
        }
        cout<<buffer<<"  ";
    }
    buffer=align_query.qual[align_query.len-1];
    if (align_target.ell[align_target.len-1]>(int)buffer.length()){
        for (int j=align_target.ell[align_target.len-1]-buffer.length(); j>0; j--){
            buffer+=" ";
        }
    }
    cout<<buffer<<endl;

    // print out target sequence
    for (int i=0; i<align_target.len-1; i++){
        buffer="";
        for (int j=0; j<align_target.ell[i]; j++){
            buffer+=align_target.alpha[i];
        }
        if (align_target.ell[i]<(int)align_query.seq[i].length()){
            for (int j=align_query.seq[i].length()-align_target.ell[i]; j>0; j--){
                buffer+=" ";
            }
        }
        cout<<buffer<<"  ";
    }
    buffer="";
    for (int j=0; j<align_target.ell[align_target.len-1]; j++){
        buffer+=align_target.alpha[align_target.len-1];
    }
    if (align_target.ell[align_target.len-1]<(int)align_query.seq[align_query.len-1].length()){
        for (int j=align_query.seq[align_query.len-1].length()-align_target.ell[align_target.len-1]; j>0; j--){
            buffer+=" ";
        }
    }
    cout<<buffer<<endl;
}


/**
 * @brief SemiHomopolymerAlignmentPool::SemiHomopolymerAlignmentPool
 */
SemiHomopolymerAlignmentPool::SemiHomopolymerAlignmentPool(){
    ;
}

/**
 * @brief SemiHomopolymerAlignmentPool::SemiHomopolymerAlignmentPool
 * @param another_pool
 */
SemiHomopolymerAlignmentPool::SemiHomopolymerAlignmentPool(const SemiHomopolymerAlignmentPool &another){
    this->align_pool=another.align_pool;
}

/**
 * @brief SemiHomopolymerAlignmentPool::~SemiHomopolymerAlignmentPool
 */
SemiHomopolymerAlignmentPool::~SemiHomopolymerAlignmentPool(){
    ;
}

/**
 * @brief SemiHomopolymerAlignmentPool::operator =
 * @param another
 * @return
 */
SemiHomopolymerAlignmentPool& SemiHomopolymerAlignmentPool::operator =(const SemiHomopolymerAlignmentPool& another){
    this->align_pool=another.align_pool;
    return *this;
}

/**
 * @brief SemiHomopolymerAlignmentPool::open
 * @param filename
 */
void SemiHomopolymerAlignmentPool::open(string filename){
    // clear everythin
    this->align_pool=vector<SemiHomopolymerAlignment>();

    // first, load in nucleotide space
    NucleotideAlignmentPool nucl_space;
    nucl_space.open(filename);
    // second, covert to homopolymer space
    for (int i=0; i<(int)nucl_space.align_pool.size(); i++){
        SemiHomopolymerAlignment aln1(nucl_space.align_pool[i]);
        this->align_pool.push_back(aln1);
    }
}

/**
 * @brief SemiHomopolymerAlignmentPool::print
 */
void SemiHomopolymerAlignmentPool::print(){
    for (int i=0; i<(int)this->align_pool.size(); i++){
        align_pool[i].print();
    }
}
