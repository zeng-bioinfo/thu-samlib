#include "nucl_align.h"

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

NucleotideAlignmentStat::~NucleotideAlignmentStat(){
    ;
}

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


//////////////////////////////////////////////////////
//      Class NucleotideAlignment                   //
//////////////////////////////////////////////////////


NucleotideAlignment::NucleotideAlignment()
{
    ;
}

NucleotideAlignment::NucleotideAlignment(const NucleotideAlignment &_align){
    this->align_name=_align.align_name;
    this->query_seq=_align.query_seq;
    this->query_qual=_align.query_qual;
    this->target_seq=_align.target_seq;
    this->align_status=_align.align_status;
}

NucleotideAlignment::~NucleotideAlignment(){
    ;
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
 * @brief NucleotideAlignmentPool::open
 * @param filename
 */
void NucleotideAlignmentPool::open(string filename){
    string tmp_name="";
    string tmp_query="";
    string tmp_qual="";
    string tmp_target="";
    string tmp_status="";

    int lc=0;
    ifstream fin(filename.c_str());
    // loop over file
    while(fin.is_open() && fin.good()){
        if (lc%5==0){
            getline(fin,tmp_name);
        }
        if (lc%5==1){
            getline(fin,tmp_query);
        }
        if (lc%5==2){
            getline(fin,tmp_qual);
        }
        if (lc%5==3){
            getline(fin,tmp_target);
        }
        if (lc%5==4){
            getline(fin,tmp_status);

            NucleotideAlignment tmp_aln;
            tmp_aln.align_name=tmp_name;
            tmp_aln.query_seq=tmp_query;
            tmp_aln.query_qual=tmp_qual;
            tmp_aln.target_seq=tmp_target;
            tmp_aln.align_status=tmp_status;
            align_pool.push_back(tmp_aln);

            tmp_name="";
            tmp_query="";
            tmp_qual="";
            tmp_target="";
            tmp_status="";
        }
        lc++;
    }
    fin.close();
}
