#include "nucl_align.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <boost/assign/list_of.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;

// TODO #1: consider the alignment status?
namespace NucleotideSpace{
    std::map<char, int> nucl2idx=boost::assign::map_list_of
            ('A',int(A))('C',int(C))('G',int(G))('T',int(T))('I',int(I))
            ('U',int(U))('W',int(W))('S',int(S))('M',int(M))('K',int(K))
            ('R',int(R))('Y',int(Y))('B',int(B))('D',int(D))('H',int(H))
            ('V',int(V))('N',int(N));
    std::map<int, char> idx2nucl=boost::assign::map_list_of
            (int(A),'A')(int(C),'C')(int(G),'G')(int(T),'T')(int(I),'I')
            (int(U),'U')(int(W),'W')(int(S),'S')(int(M),'M')(int(K),'K')
            (int(R),'R')(int(Y),'Y')(int(B),'B')(int(D),'D')(int(H),'H')
            (int(V),'V')(int(N),'N');
}

//int NucleotideSpace::index(char c){
//    using namespace NucleotideSpace;
//    int nucl;

//    switch (c){
//    case 'A':
//        nucl=A;
//        break;
//    case 'C':
//        nucl=C;
//        break;
//    case 'G':
//        nucl=G;
//        break;
//    case 'T':
//        nucl=T;
//        break;
//    case '-':
//        nucl=I;
//        break;
//    case 'U':
//        nucl=U;
//        break;
//    case 'W':
//        nucl=W;
//        break;
//    case 'S':
//        nucl=S;
//        break;
//    case 'M':
//        nucl=M;
//        break;
//    case 'K':
//        nucl=K;
//        break;
//    case 'R':
//        nucl=R;
//        break;
//    case 'Y':
//        nucl=Y;
//        break;
//    case 'B':
//        nucl=B;
//        break;
//    case 'D':
//        nucl=D;
//        break;
//    case 'H':
//        nucl=H;
//        break;
//    case 'V':
//        nucl=V;
//        break;
//    case 'N':
//        nucl=N;
//        break;
//    }

//    return nucl;
//}

//string NucleotideSpace::symbol(Nucleotide nucl){
//    string letter;
//    switch(nucl){
//    case A:
//        letter="A";
//        break;
//    case C:
//        letter="C";
//        break;
//    case G:
//        letter="G";
//        break;
//    case T:
//        letter="T";
//        break;
//    case I:
//        letter="-";
//        break;
//    }
//    return letter;
//}

//////////////////////////////////////////////////////
//      Class NucleotideAlignment Statistics        //
//////////////////////////////////////////////////////

/**
 * @brief NucleotideAlignmentStat::NucleotideAlignmentStat
 */
NucleotideAlignmentStat::NucleotideAlignmentStat(){
    // 0-order statistics
    // set base_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]=0;
        }
    }
    // set qual_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]=0;
        }
    }
    // set base_qual_cooc_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::QualityScoreMax; k++){
                this->base_qual_cooc_table[i][j][k]=0;
            }
        }
    }

    // set gc_base_qual_count_table to 0
    for (int i=0; i<100; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::NucleotideSize; k++){
                for (int q=0; q<NucleotideSpace::QualityScoreMax; q++){
                    this->gc_base_qual_count_table[i][j][k][q]=0;
                }
            }
        }
    }
}

/**
 * @brief NucleotideAlignmentStat::NucleotideAlignmentStat
 * @param _align_stat
 */
NucleotideAlignmentStat::NucleotideAlignmentStat(const NucleotideAlignmentStat &_align_stat){
    // 0-order statistics
    // set base_call_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]=_align_stat.base_call_table[i][j];
        }
    }
    // set qual_call_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]=_align_stat.qual_call_table[i][j];
        }
    }
    // set base_qual_cooc_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::QualityScoreMax; k++){
                this->base_qual_cooc_table[i][j][k]=_align_stat.base_qual_cooc_table[i][j][k];
            }
        }
    }
    // set gc_base_qual_count_table to given item
    for (int i=0; i<100; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::NucleotideSize; k++){
                for (int q=0; q<NucleotideSpace::QualityScoreMax; q++){
                    this->gc_base_qual_count_table[i][j][k][q]=_align_stat.gc_base_qual_count_table[i][j][k][q];
                }
            }
        }
    }
}

/**
 * @brief NucleotideAlignmentStat::~NucleotideAlignmentStat
 */
NucleotideAlignmentStat::~NucleotideAlignmentStat(){
    ;
}

/**
 * @brief NucleotideAlignmentStat::operator =
 * @param _align_stat
 * @return
 */
NucleotideAlignmentStat& NucleotideAlignmentStat::operator =(const NucleotideAlignmentStat& _align_stat){
    // 0-order statistics
    // set base_call_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]=_align_stat.base_call_table[i][j];
        }
    }
    // set qual_call_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]=_align_stat.qual_call_table[i][j];
        }
    }
    // set base_qual_cooc_table to given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::QualityScoreMax; k++){
                this->base_qual_cooc_table[i][j][k]=_align_stat.base_qual_cooc_table[i][j][k];
            }
        }
    }
    // set gc_base_qual_count_table to given item
    for (int i=0; i<100; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::NucleotideSize; k++){
                for (int q=0; q<NucleotideSpace::QualityScoreMax; q++){
                    this->gc_base_qual_count_table[i][j][k][q]=_align_stat.gc_base_qual_count_table[i][j][k][q];
                }
            }
        }
    }

    return *this;
}

/**
 * @brief NucleotideAlignmentStat::operator +=
 * @param _align_stat
 * @return
 */
NucleotideAlignmentStat& NucleotideAlignmentStat::operator +=(const NucleotideAlignmentStat& _align_stat){
    // 0-order statistics
    // update base_call_table by given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]+=_align_stat.base_call_table[i][j];
        }
    }
    // update qual_call_table by given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]+=_align_stat.qual_call_table[i][j];
        }
    }
    // update base_qual_cooc_table by given item
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::QualityScoreMax; k++){
                this->base_qual_cooc_table[i][j][k]+=_align_stat.base_qual_cooc_table[i][j][k];
            }
        }
    }
    // update gc_base_qual_count_table by given item
    for (int i=0; i<100; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            for (int k=0; k<NucleotideSpace::NucleotideSize; k++){
                for (int q=0; q<NucleotideSpace::QualityScoreMax; q++){
                    this->gc_base_qual_count_table[i][j][k][q]+=_align_stat.gc_base_qual_count_table[i][j][k][q];
                }
            }
        }
    }

    return *this;
}

/**
 * @brief NucleotideAlignmentStat::print
 * @param c
 */
void NucleotideAlignmentStat::print(int c){
    using namespace NucleotideSpace;
    char buffer[8192];

    if (c==0){
        // header
        sprintf(buffer, "cycle gc qual match mismatch error");
        for (int alpha=(int)A; alpha<=(int)T; alpha++){
            for (int beta=(int)A; beta<=(int)T; beta++){
                sprintf(buffer,"%s %c%c", buffer, idx2nucl[alpha], idx2nucl[beta]);
            }
        }
        sprintf(buffer, "%s %s", buffer, "insert");
        for (int beta=(int)A; beta<=(int)T; beta++){
            sprintf(buffer, "%s %c%c", buffer, '-', idx2nucl[beta]);
        }
        sprintf(buffer, "%s %s", buffer, "delete");
        for (int alpha=(int)A; alpha<=(int)T; alpha++){
            sprintf(buffer, "%s %c%c", buffer, idx2nucl[alpha], '-');
        }
        cout<<buffer<<endl;
    }

    // print out
    for (int gc=0; gc<100; gc+=10){
        for (int qual=0; qual<100; qual+=10){
            sprintf(buffer, "%d %d %d", c+1, gc/10, qual/10);
            // match, mismatch, and error rate
            int match=0;
            int mismatch=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int a=(int)A; a<=(int)T; a++){
                        for (int b=(int)A; b<=(int)T; b++){
                            if (a!=b){
                                mismatch+=gc_base_qual_count_table[k][a][b][q];
                            }else{
                                match+=gc_base_qual_count_table[k][a][b][q];
                            }
                        }
                    }
                }
            }
            sprintf(buffer, "%s %d %d %f", buffer, match, mismatch, mismatch/(match+1e-9));
            // a->b
            for (int a=(int)A; a<=(int)T; a++){
                for (int b=(int)A; b<=(int)T; b++){
                    int count=0;
                    for (int k=gc; k<gc+10; k++){
                        for (int q=qual; q<qual+10; q++){
                            count+=gc_base_qual_count_table[k][a][b][q];
                        }
                    }
                    sprintf(buffer, "%s %d", buffer, count);
                }
            }
            // insertion
            int ins=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int b=(int)A; b<=(int)T; b++){
                        ins+=gc_base_qual_count_table[k][(int)I][b][q];
                    }
                }
            }
            sprintf(buffer, "%s %d", buffer, ins);
            // \-->b
            for (int b=(int)A; b<=(int)T; b++){
                int count=0;
                for (int k=gc; k<gc+10; k++){
                    for (int q=qual; q<qual+10; q++){
                        count+=gc_base_qual_count_table[k][(int)I][b][q];
                    }
                }
                sprintf(buffer, "%s %d", buffer, count);
            }
            // deletion
            int del=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int a=(int)A; a<=(int)T; a++){
                        del+=gc_base_qual_count_table[k][a][(int)I][q];
                    }
                }
            }
            sprintf(buffer, "%s %d", buffer, del);
            // a->\-
            for (int a=(int)A; a<=(int)T; a++){
                int count=0;
                for (int k=gc; k<gc+10; k++){
                    for (int q=qual; q<qual+10; q++){
                        count+=gc_base_qual_count_table[k][a][(int)I][q];
                    }
                }
                sprintf(buffer, "%s %d", buffer, count);
            }
            cout<<buffer<<endl;
        }
    }
}

/**
 * @brief NucleotideAlignmentStat::print
 * @param ofs
 */
void NucleotideAlignmentStat::print(ofstream& ofs, int c){
    // TODO: output to the specified file
    using namespace NucleotideSpace;
    char buffer[8192];

    if (c==0){
        // header
        sprintf(buffer, "cycle gc qual match mismatch error");
        for (int alpha=(int)A; alpha<=(int)T; alpha++){
            for (int beta=(int)A; beta<=(int)T; beta++){
                sprintf(buffer,"%s %c%c", buffer, idx2nucl[alpha], idx2nucl[beta]);
            }
        }
        sprintf(buffer, "%s %s", buffer, "insert");
        for (int beta=(int)A; beta<=(int)T; beta++){
            sprintf(buffer, "%s %c%c", buffer, '-', idx2nucl[beta]);
        }
        sprintf(buffer, "%s %s", buffer, "delete");
        for (int alpha=(int)A; alpha<=(int)T; alpha++){
            sprintf(buffer, "%s %c%c", buffer, idx2nucl[alpha], '-');
        }
        ofs<<buffer<<endl;
    }

    // print out
    for (int gc=0; gc<100; gc+=10){
        for (int qual=0; qual<100; qual+=10){
            sprintf(buffer, "%d %d %d", c+1, gc/10, qual/10);
            // match, mismatch, and error rate
            int match=0;
            int mismatch=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int a=(int)A; a<=(int)T; a++){
                        for (int b=(int)A; b<=(int)T; b++){
                            if (a!=b){
                                mismatch+=gc_base_qual_count_table[k][a][b][q];
                            }else{
                                match+=gc_base_qual_count_table[k][a][b][q];
                            }
                        }
                    }
                }
            }
            sprintf(buffer, "%s %d %d %f", buffer, match, mismatch, mismatch/(match+1e-9));
            // a->b
            for (int a=(int)A; a<=(int)T; a++){
                for (int b=(int)A; b<=(int)T; b++){
                    int count=0;
                    for (int k=gc; k<gc+10; k++){
                        for (int q=qual; q<qual+10; q++){
                            count+=gc_base_qual_count_table[k][a][b][q];
                        }
                    }
                    sprintf(buffer, "%s %d", buffer, count);
                }
            }
            // insertion
            int ins=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int b=(int)A; b<=(int)T; b++){
                        ins+=gc_base_qual_count_table[k][(int)I][b][q];
                    }
                }
            }
            sprintf(buffer, "%s %d", buffer, ins);
            // \-->b
            for (int b=(int)A; b<=(int)T; b++){
                int count=0;
                for (int k=gc; k<gc+10; k++){
                    for (int q=qual; q<qual+10; q++){
                        count+=gc_base_qual_count_table[k][(int)I][b][q];
                    }
                }
                sprintf(buffer, "%s %d", buffer, count);
            }
            // deletion
            int del=0;
            for (int k=gc; k<gc+10; k++){
                for (int q=qual; q<qual+10; q++){
                    for (int a=(int)A; a<=(int)T; a++){
                        del+=gc_base_qual_count_table[k][a][(int)I][q];
                    }
                }
            }
            sprintf(buffer, "%s %d", buffer, del);
            // a->\-
            for (int a=(int)A; a<=(int)T; a++){
                int count=0;
                for (int k=gc; k<gc+10; k++){
                    for (int q=qual; q<qual+10; q++){
                        count+=gc_base_qual_count_table[k][a][(int)I][q];
                    }
                }
                sprintf(buffer, "%s %d", buffer, count);
            }
            ofs<<buffer<<endl;
        }
    }
}

//////////////////////////////////////////////////////
//      Class NucleotideAlignment                   //
//////////////////////////////////////////////////////

/**
 * @brief NucleotideAlignment::NucleotideAlignment
 */
NucleotideAlignment::NucleotideAlignment()
{
    this->align_name="";
    this->query_seq="";
    this->query_qual="";
    this->target_seq="";
    this->align_status="";
}

/**
 * @brief NucleotideAlignment::NucleotideAlignment
 * @param _align
 */
NucleotideAlignment::NucleotideAlignment(const NucleotideAlignment &_align){
    this->align_name=_align.align_name;
    this->query_seq=_align.query_seq;
    this->query_qual=_align.query_qual;
    this->query_strand=_align.query_strand;
    this->target_seq=_align.target_seq;
    this->align_status=_align.align_status;

    this->raw_query=_align.raw_query;
    this->raw_quality=_align.raw_quality;
    this->raw_target=_align.raw_target;
}

/**
 * @brief NucleotideAlignment::~NucleotideAlignment
 */
NucleotideAlignment::~NucleotideAlignment(){
    ;
}

/**
 * @brief NucleotideAlignment::operator =
 * @param _align
 * @return
 */
NucleotideAlignment& NucleotideAlignment::operator =(const NucleotideAlignment& _align){
    this->align_name=_align.align_name;
    this->query_seq=_align.query_seq;
    this->query_qual=_align.query_qual;
    this->query_strand=_align.query_strand;
    this->target_seq=_align.target_seq;
    this->align_status=_align.align_status;

    this->raw_query=_align.raw_query;
    this->raw_quality=_align.raw_quality;
    this->raw_target=_align.raw_target;

    return *this;
}

/**
 * @brief NucleotideAlignment::indel_shift_right
 */
void NucleotideAlignment::indel_shift_right(){
    // scan over the alignment
    for (int i=1; i<(int)query_seq.length(); i++){       // jump over
        if (query_seq[i]=='-' || target_seq[i]=='-' || query_seq[i]!=target_seq[i]){
            continue;
        }else{                                      // query_seq[i]==target_seq[i]
            for (int j=i-1; j>=0; j--){
                if (query_seq[j]!='-' && target_seq[j]!='-'){
                    break;
                }else{
                    if (query_seq[j]=='-' && target_seq[j]==query_seq[j+1]){
                        // swap query_seq[j] and query_seq[j+1]
                        std::swap(query_seq[j], query_seq[j+1]);
                        // swap query_qual[j] and query_qual[j+1]
                        std::swap(query_qual[j], query_qual[j+1]);
                        // swap align_status[j] and align_status[j+1]
                        std::swap(align_status[j], align_status[j+1]);
                    }
                    else if (target_seq[j]=='-' && query_seq[j]==target_seq[j+1]){
                        // swap target_seq[j] and target_seq[j+1]
                        std::swap(target_seq[j], target_seq[j+1]);
                        // swap align_status[j] and align_status[j+1]
                        std::swap(align_status[j], align_status[j+1]);
                    }
                }
            }
        }
    }
}


void NucleotideAlignment::statistics(vector<NucleotideAlignmentStat> &stat, int cycles){
    using namespace NucleotideSpace;
    // initialze or re-size the vector
    if (!stat.empty()) {
        stat.clear();
    }
    for (int c=0; c<cycles; c++){
        NucleotideAlignmentStat s;
        stat.push_back(s);
    }

    // compute the length of the target sequence
    int len=0;
    for (int i=0; i<(int)target_seq.length(); i++){
        if (target_seq[i]=='-') {
            continue;
        }
        len++;
    }

    // compute gc percentage
    int gc=0;
    int atgc=0;
    for (int i=0; i<query_seq.length(); i++){
        if (query_seq[i]=='-') continue;

        if (query_seq[i]=='G' || query_seq[i]=='C')
            gc+=1;
        atgc+=1;
    }
    gc=int(gc*100/(atgc+0.));

    // scan the alignment to compute the statistics
    for (int i=0,j=0; i<(int)align_status.length(); i++){    // i is a pointer on the alignment, j is a pointer on the target
        // computer the cycle no
        int c=(int)(j/(len+0.)*cycles);
        // statistics object
        NucleotideAlignmentStat s=stat[c];
        // collect the statistics
        if (nucl2idx[target_seq[i]]<NucleotideSize && nucl2idx[query_seq[i]]<NucleotideSize){
            // base_call count table
            s.base_call_table[nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]]+=1;
            // qual_call_count table
            if (query_seq[i]!='-'){
                s.qual_call_table[nucl2idx[target_seq[i]]][(int)query_qual[i]-33]+=1;
                s.base_qual_cooc_table[nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]][(int)query_qual[i]-33]+=1;
            }else{
                s.qual_call_table[nucl2idx[target_seq[i]]][0]+=1;
            }
            // base_qual_count table
            if (query_seq[i]!='-'){
                s.base_qual_cooc_table[nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]][(int)query_qual[i]-33]+=1;
            }else{
                s.base_qual_cooc_table[nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]][0]+=1;
            }
            // gc_base_qual_count table
            if (query_seq[i]!='-'){
                s.gc_base_qual_count_table[gc][nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]][(int)query_qual[i]-33]+=1;
            }else{
                s.gc_base_qual_count_table[gc][nucl2idx[target_seq[i]]][nucl2idx[query_seq[i]]][0]+=1;
            }
        }
        // restore the statistics object
        stat[c]=s;
        // update j
        if (target_seq[i]!='-'){
            j+=1;
        }
    }
}

//////////////////////////////////////////////////////
//      Class NucleotideAlignment  Pool             //
//////////////////////////////////////////////////////

/**
 * @brief NucleotideAlignmentPool::NucleotideAlignmentPool
 */
NucleotideAlignmentPool::NucleotideAlignmentPool(){
    ;
}

/**
 * @brief NucleotideAlignmentPool::~NucleotideAlignmentPool
 */
NucleotideAlignmentPool::~NucleotideAlignmentPool(){
    ;
}

/**
 * @brief NucleotideAlignmentPool::statistics
 * @param cycles
 */
void NucleotideAlignmentPool::statistics(int cycles){
    vector<NucleotideAlignmentStat> stat;
    statistics(stat,cycles);
    // print out
    for (int c=0; c<cycles; c++){
        NucleotideAlignmentStat s=stat[c];
        s.print(c);
    }
//    using namespace NucleotideSpace;
//    char buffer[8192];
//    // 0-order statistics
//    // title
//    sprintf(buffer, "cycle gc qual match mismatch error");
//    for (int alpha=(int)A; alpha<=(int)T; alpha++){
//        for (int beta=(int)A; beta<=(int)T; beta++){
//            sprintf(buffer,"%s %c%c", buffer, idx2nucl[alpha], idx2nucl[beta]);
//        }
//    }
//    sprintf(buffer, "%s %s", buffer, "insert");
//    for (int beta=(int)A; beta<=(int)T; beta++){
//        sprintf(buffer, "%s %c%c", buffer, '-', idx2nucl[beta]);
//    }
//    sprintf(buffer, "%s %s", buffer, "delete");
//    for (int alpha=(int)A; alpha<=(int)T; alpha++){
//        sprintf(buffer, "%s %c%c", buffer, idx2nucl[alpha], '-');
//    }
//    cout<<buffer<<endl;
//    // print out statistics
//    for (int c=0; c<cycles; c++){
//        NucleotideAlignmentStat s=stat[c];
//        for (int gc=0; gc<100; gc+=10){
//            for (int qual=0; qual<100; qual+=10){
//                sprintf(buffer, "%d %d %d", c+1, gc/10, qual/10);
//                // match, mismatch, and error rate
//                int match=0;
//                int mismatch=0;
//                for (int k=gc; k<gc+10; k++){
//                    for (int q=qual; q<qual+10; q++){
//                        for (int a=(int)A; a<=(int)T; a++){
//                            for (int b=(int)A; b<=(int)T; b++){
//                                if (a!=b){
//                                    mismatch+=s.gc_base_qual_count_table[k][a][b][q];
//                                }else{
//                                    match+=s.gc_base_qual_count_table[k][a][b][q];
//                                }
//                            }
//                        }
//                    }
//                }
//                sprintf(buffer, "%s %d %d %f", buffer, match, mismatch, mismatch/(match+1e-9));
//                // a->b
//                for (int a=(int)A; a<=(int)T; a++){
//                    for (int b=(int)A; b<=(int)T; b++){
//                        int count=0;
//                        for (int k=gc; k<gc+10; k++){
//                            for (int q=qual; q<qual+10; q++){
//                                count+=s.gc_base_qual_count_table[k][a][b][q];
//                            }
//                        }
//                        sprintf(buffer, "%s %d", buffer, count);
//                    }
//                }
//                // insertion
//                int ins=0;
//                for (int k=gc; k<gc+10; k++){
//                    for (int q=qual; q<qual+10; q++){
//                        for (int b=(int)A; b<=(int)T; b++){
//                            ins+=s.gc_base_qual_count_table[k][(int)I][b][q];
//                        }
//                    }
//                }
//                sprintf(buffer, "%s %d", buffer, ins);
//                // \-->b
//                for (int b=(int)A; b<=(int)T; b++){
//                    int count=0;
//                    for (int k=gc; k<gc+10; k++){
//                        for (int q=qual; q<qual+10; q++){
//                            count+=s.gc_base_qual_count_table[k][(int)I][b][q];
//                        }
//                    }
//                    sprintf(buffer, "%s %d", buffer, count);
//                }
//                // deletion
//                int del=0;
//                for (int k=gc; k<gc+10; k++){
//                    for (int q=qual; q<qual+10; q++){
//                        for (int a=(int)A; a<=(int)T; a++){
//                            del+=s.gc_base_qual_count_table[k][a][(int)I][q];
//                        }
//                    }
//                }
//                sprintf(buffer, "%s %d", buffer, del);
//                // a->\-
//                for (int a=(int)A; a<=(int)T; a++){
//                    int count=0;
//                    for (int k=gc; k<gc+10; k++){
//                        for (int q=qual; q<qual+10; q++){
//                            count+=s.gc_base_qual_count_table[k][a][(int)I][q];
//                        }
//                    }
//                    sprintf(buffer, "%s %d", buffer, count);
//                }
//                cout<<buffer<<endl;
//            }
//        }
//    }
}

/**
 * @brief NucleotideAlignmentPool::statistics
 * @param outfile
 * @param cycles
 */
void NucleotideAlignmentPool::statistics(string outfile, int cycles){
    vector<NucleotideAlignmentStat> stat;
    statistics(stat,cycles);

    ofstream ofs;
    ofs.open(outfile);
    // print out
    for (int c=0; c<cycles; c++){
        NucleotideAlignmentStat s=stat[c];
        s.print(ofs, c);
    }
    ofs.close();
}

/**
 * @brief NucleotideAlignmentPool::statistics
 * @param stat
 * @param cycles
 */
void NucleotideAlignmentPool::statistics(vector<NucleotideAlignmentStat> &stat, int cycles){
    // initialize stat
    if (!stat.empty()){
        stat.clear();
    }
    for (int c=0; c<cycles; c++){
        NucleotideAlignmentStat s;
        stat.push_back(s);
    }

    // loop over all alignments
    for (int i=0; i<(int)align_pool.size(); i++){
        NucleotideAlignment aln1=align_pool[i];
        // the statistics of the current alignment
        vector<NucleotideAlignmentStat> aln1_stat;
        aln1.statistics(aln1_stat,cycles);
        // update global statistics
        for (int c=0; c<(int)stat.size(); c++){
            stat[c]+=aln1_stat[c];
        }
    }
}

/**
 * @brief NucleotideAlignmentPool::open
 * @param filename
 */
void NucleotideAlignmentPool::open(string filename){
    string tmp_name="";
    string tmp_query="";
    string tmp_qual="";
    int tmp_strand=0;
    string tmp_target="";
    string tmp_status="";
    string tmp_raw_query="";
    string tmp_raw_quality="";
    string tmp_raw_target="";

    int lc=0;
    ifstream fin(filename.c_str());
    // loop over file
    while(fin.is_open() && fin.good()){
        if (lc%5==0){
            getline(fin,tmp_name);

            // strip newline character
            if (*tmp_name.rbegin()=='\r'){
                tmp_name.erase(tmp_name.length()-1);
            }

            // extract the strand information
            size_t pos =tmp_name.find("strand:");
            if (pos!=string::npos){
                string strand=tmp_name.substr(pos+7,1);
                if (strand=="+"){
                    tmp_strand=1;
                }else{
                    tmp_strand=0;
                }
            }
        }
        if (lc%5==1){
            getline(fin,tmp_query);

            // strip newline character
            if (*tmp_query.rbegin()=='\r'){
                tmp_query.erase(tmp_query.length()-1);
            }
        }
        if (lc%5==2){
            getline(fin,tmp_qual);

            // strip newline character
            if (*tmp_qual.rbegin()=='\r'){
                tmp_qual.erase(tmp_qual.length()-1);
            }
        }
        if (lc%5==3){
            getline(fin,tmp_target);
            // strip newline character
            if (*tmp_target.rbegin()=='\r'){
                tmp_target.erase(tmp_target.length()-1);
            }
        }
        if (lc%5==4){
            getline(fin,tmp_status);

            // strip newline character
            if (*tmp_status.rbegin()=='\r'){
                tmp_status.erase(tmp_status.length()-1);
            }

            // raw query and quality
            tmp_raw_query="";
            tmp_raw_quality="";
            for (int t=0; t<(int)tmp_query.length(); t++){
                if (tmp_query[t]!='-' && tmp_query[t]!='\r'){
                    tmp_raw_query+=tmp_query[t];
                    tmp_raw_quality+=tmp_qual[t];
                }
            }
            // raw target
            tmp_raw_target="";
            for (int t=0; t<(int)tmp_target.length(); t++){
                if (tmp_target[t]!='-' && tmp_target[t]!='\r' &&
                        tmp_status[t]!='S' && tmp_status[t]!='s'){
                    tmp_raw_target+=tmp_target[t];
                }
            }

            // package up an alignment
            NucleotideAlignment tmp_aln;
            tmp_aln.align_name=tmp_name;
            tmp_aln.query_seq=tmp_query;
            tmp_aln.query_qual=tmp_qual;
            tmp_aln.query_strand=tmp_strand;
            tmp_aln.target_seq=tmp_target;
            tmp_aln.align_status=tmp_status;
            tmp_aln.raw_query=tmp_raw_query;
            tmp_aln.raw_quality=tmp_raw_quality;
            tmp_aln.raw_target=tmp_raw_target;
            tmp_aln.indel_shift_right();
            align_pool.push_back(tmp_aln);

            // clear everything
            tmp_name="";
            tmp_query="";
            tmp_qual="";
            tmp_strand=0;
            tmp_target="";
            tmp_status="";
            tmp_raw_query="";
            tmp_raw_quality="";
            tmp_raw_target="";

//            // debug
//            int t0,t1,q0,q1;
//            NucleotideAlignmentMethod alignment_band;
//            alignment_band.find_exact_match_segment_chain(tmp_aln,
//                    t0,t1,q0,q1);


//            cout<<tmp_aln.align_name<<endl;
//            cout<<t0<<" <-> "<<q0<<endl
//                <<t1<<" <-> "<<q1<<endl;
//            cout<<endl;
        }
        lc++;
    }
    fin.close();
}
