#include "nucl_align.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <algorithm>
using namespace std;

// TODO #1: consider the alignment status?

int NucleotideSpace::index(char c){
    using namespace NucleotideSpace;
    int nucl;

    switch (c){
    case 'A':
        nucl=A;
        break;
    case 'C':
        nucl=C;
        break;
    case 'G':
        nucl=G;
        break;
    case 'T':
        nucl=T;
        break;
    case '-':
        nucl=I;
        break;
    case 'U':
        nucl=U;
        break;
    case 'W':
        nucl=W;
        break;
    case 'S':
        nucl=S;
        break;
    case 'M':
        nucl=M;
        break;
    case 'K':
        nucl=K;
        break;
    case 'R':
        nucl=R;
        break;
    case 'Y':
        nucl=Y;
        break;
    case 'B':
        nucl=B;
        break;
    case 'D':
        nucl=D;
        break;
    case 'H':
        nucl=H;
        break;
    case 'V':
        nucl=V;
        break;
    case 'N':
        nucl=N;
        break;
    }

    return nucl;
}

string NucleotideSpace::symbol(Nucleotide nucl){
    string letter;
    switch(nucl){
    case A:
        letter="A";
        break;
    case C:
        letter="C";
        break;
    case G:
        letter="G";
        break;
    case T:
        letter="T";
        break;
    case I:
        letter="-";
        break;
    }
    return letter;
}

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
}

/**
 * @brief NucleotideAlignmentStat::NucleotideAlignmentStat
 * @param _align_stat
 */
NucleotideAlignmentStat::NucleotideAlignmentStat(const NucleotideAlignmentStat &_align_stat){
    // 0-order statistics
    // set base_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]=_align_stat.base_call_table[i][j];
        }
    }
    // set qual_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]=_align_stat.qual_call_table[i][j];
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
    // set base_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]=_align_stat.base_call_table[i][j];
        }
    }
    // set qual_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]=_align_stat.qual_call_table[i][j];
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
    // set base_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::NucleotideSize; j++){
            this->base_call_table[i][j]+=_align_stat.base_call_table[i][j];
        }
    }
    // set qual_call_table to 0
    for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
        for (int j=0; j<NucleotideSpace::QualityScoreMax; j++){
            this->qual_call_table[i][j]+=_align_stat.qual_call_table[i][j];
        }
    }

    return *this;
}

/**
 * @brief NucleotideAlignmentStat::print
 */
void NucleotideAlignmentStat::print(){
    using namespace NucleotideSpace;
    // 0-order statistics
    // base call table
    cout<<"[Base-Call Table]"<<endl;
    printf("Call-To-Be  %10s %10s %10s %10s %10s\n",
           symbol(A).c_str(),symbol(C).c_str(),symbol(G).c_str(),symbol(T).c_str(),symbol(I).c_str());
    for (int nucl=A; nucl<NucleotideSize; nucl++){
        printf("Input    %s  %10d %10d %10d %10d %10d\n",symbol((Nucleotide)nucl).c_str(),this->base_call_table[nucl][A],
               this->base_call_table[nucl][C],this->base_call_table[nucl][G],this->base_call_table[nucl][T],
               this->base_call_table[nucl][I]);
    }

    // quality call table
    char buffer[8192];
    cout<<"[Quality-Score-Call Table]"<<endl;
    sprintf(buffer, "Quality-Score ");
    for (int nucl=A; nucl<NucleotideSize; nucl++){
        sprintf(buffer,"%s %10s",buffer,symbol((Nucleotide)nucl).c_str());
    }
    cout<<buffer<<endl;
    for (int q=0; q<QualityScoreMax; q++){
        sprintf(buffer,"%2d            ",q);
        for (int nucl=A; nucl<NucleotideSize; nucl++){
            sprintf(buffer,"%s %10d",buffer,this->qual_call_table[nucl][q]);
        }
        cout<<buffer<<endl;
    }
}

/**
 * @brief NucleotideAlignmentStat::print
 * @param ofs
 */
void NucleotideAlignmentStat::print(ofstream& ofs){
    using namespace NucleotideSpace;
    char buffer[8192];
    // 0-order statistics
    // base call table
    ofs<<"[Base-Call Table]"<<endl;
    sprintf(buffer,"Call-To-Be ");
    for (int nucl=A; nucl<NucleotideSize; nucl++){
        sprintf(buffer,"%s %10s", buffer, symbol((Nucleotide)nucl).c_str());
    }
    ofs<<buffer<<endl;
    for (int nucl=A; nucl<NucleotideSize; nucl++){
        sprintf(buffer,"Input    %s ",symbol((Nucleotide)nucl).c_str());
        for (int nucl2=A; nucl2<NucleotideSize; nucl2++){
            sprintf(buffer,"%s %10d", buffer, this->base_call_table[nucl][nucl2]);
        }
        ofs<<buffer<<endl;
    }
    // quality call table
    ofs<<"[Quality-Score-Call Table]"<<endl;
    sprintf(buffer, "Quality-Score ");
    for (int nucl=A; nucl<NucleotideSize; nucl++){
        sprintf(buffer,"%s %10s",buffer,symbol((Nucleotide)nucl).c_str());
    }
    ofs<<buffer<<endl;
    for (int q=0; q<QualityScoreMax; q++){
        sprintf(buffer,"%2d            ",q);
        for (int nucl=A; nucl<NucleotideSize; nucl++){
            sprintf(buffer,"%s %10d",buffer,this->qual_call_table[nucl][q]);
        }
        ofs<<buffer<<endl;
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

    // scan the alignment to compute the statistics
    for (int i=0,j=0; i<(int)align_status.length(); i++){    // i is a pointer on the alignment, j is a pointer on the target
        // computer the cycle no
        int c=(int)(j/(len+0.)*cycles);
        // statistics object
        NucleotideAlignmentStat s=stat[c];
        // collect the statistics
        if (index(target_seq[i])<NucleotideSize && index(query_seq[i])<NucleotideSize){
            s.base_call_table[index(target_seq[i])][index(query_seq[i])]+=1;
            s.qual_call_table[index(target_seq[i])][(int)query_qual[i]-33]+=1;
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
        cout<<"[Cycle #]"<<c<<endl;
        NucleotideAlignmentStat s=stat[c];
        s.print();
        cout<<endl<<endl;
    }
}

/**
 * @brief NucleotideAlignmentPool::statistics
 * @param outfile
 * @param cycles
 */
void NucleotideAlignmentPool::statistics(string outfile, int cycles){
    vector<NucleotideAlignmentStat> stat;
    statistics(stat,cycles);
    // print out
    for (int c=0; c<cycles; c++){
        string of=outfile+"."+to_string(c);
        ofstream ofs(of.c_str());
        NucleotideAlignmentStat s=stat[c];
        s.print(ofs);
        ofs.close();
    }
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

            for (int t=0; t<(int)tmp_query.length(); t++){
                if (tmp_query[t]!='-' && tmp_query[t]!='\r'){
                    tmp_raw_query+=tmp_query[t];
                }
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

            for (int t=0; t<(int)tmp_target.length(); t++){
                if (tmp_target[t]!='-' && tmp_target[t]!='\r'){
                    tmp_raw_target+=tmp_target[t];
                }
            }
        }
        if (lc%5==4){
            getline(fin,tmp_status);

            // strip newline character
            if (*tmp_status.rbegin()=='\r'){
                tmp_status.erase(tmp_status.length()-1);
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
            tmp_raw_target="";

            int t0,t1,q0,q1;
            NucleotideAlignmentMethod alignment_band;
            alignment_band.find_exact_match_segment_chain(tmp_aln,
                    t0,t1,q0,q1);
            cout<<tmp_aln.align_name<<endl;
            cout<<t0<<" <-> "<<q0<<endl
                <<t1<<" <-> "<<q1<<endl;
            cout<<endl;
        }
        lc++;
    }
    fin.close();
}
