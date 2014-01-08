#include "semi_homo_align.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
using namespace std;

/**
 * @brief SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat
 */
SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat(){
    // base call table
    this->base_call_table=vector<int>(1<<(HomopolymerSpace::STRANDBITSIZE+
                                          HomopolymerSpace::ALPHABITSIZE+
                                          HomopolymerSpace::LENGTHBITSIZE+
                                          HomopolymerSpace::BETABITSIZE));
    // quality score call table
    this->qual_call_table=vector<int>(1<<(HomopolymerSpace::STRANDBITSIZE+
                                          HomopolymerSpace::ALPHABITSIZE+
                                          HomopolymerSpace::LENGTHBITSIZE+
                                          HomopolymerSpace::QUALITYBITSIZE));
    // length call table
    this->len_call_table=vector<int>(1<<(HomopolymerSpace::STRANDBITSIZE+
                                         HomopolymerSpace::ALPHABITSIZE+
                                         HomopolymerSpace::LENGTHBITSIZE+
                                         HomopolymerSpace::LENGTHBITSIZE));
}

/**
 * @brief SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat
 * @param another
 */
SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat(const SemiHomopolymerAlignmentStat &another){
    this->base_call_table=another.base_call_table;
    this->qual_call_table=another.qual_call_table;
    this->len_call_table=another.len_call_table;
}

/**
 * @brief SemiHomopolymerAlignmentStat::~SemiHomopolymerAlignmentStat
 */
SemiHomopolymerAlignmentStat::~SemiHomopolymerAlignmentStat(){
    this->base_call_table=vector<int>();
    this->qual_call_table=vector<int>();
    this->len_call_table=vector<int>();
}

/**
 * @brief SemiHomopolymerAlignmentStat::operator =
 * @param another
 * @return
 */
SemiHomopolymerAlignmentStat& SemiHomopolymerAlignmentStat::operator =(const SemiHomopolymerAlignmentStat &another){
    this->base_call_table=another.base_call_table;
    this->qual_call_table=another.qual_call_table;
    this->len_call_table=another.len_call_table;
    return *this;
}

/**
 * @brief SemiHomopolymerAlignmentStat::operator +=
 * @param another
 * @return
 */
SemiHomopolymerAlignmentStat& SemiHomopolymerAlignmentStat::operator +=(const SemiHomopolymerAlignmentStat &another){
    // update base call table
    for (int i=0; i<(int)base_call_table.size(); i++){
        base_call_table[i]+=another.base_call_table[i];
    }
    // update quality score call table
    for (int i=0; i<(int)qual_call_table.size(); i++){
        qual_call_table[i]+=another.qual_call_table[i];
    }
    // update length call table
    for (int i=0; i<(int)len_call_table.size(); i++){
        len_call_table[i]+=another.len_call_table[i];
    }

    return *this;
}

/**
 * @brief SemiHomopolymerAlignmentStat::bc_incr1
 * @param alpha
 * @param ell
 * @param beta
 */
void SemiHomopolymerAlignmentStat::bc_incr1(int strand, char alpha, int ell, char beta){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE)) |
                 (ell<<HomopolymerSpace::BETABITSIZE) | (HomopolymerSpace::index(beta));
    base_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::qc_incr1
 * @param alpha
 * @param ell
 * @param qual
 */
void SemiHomopolymerAlignmentStat::qc_incr1(int strand, char alpha, int ell, int qual){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE)) |
                 (ell<<HomopolymerSpace::QUALITYBITSIZE) | qual;
    qual_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::lc_incr1
 * @param alpha
 * @param ell
 * @param kappa
 */
void SemiHomopolymerAlignmentStat::lc_incr1(int strand, char alpha, int ell, int kappa){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE)) |
                 (ell<<HomopolymerSpace::LENGTHBITSIZE) | kappa;
    len_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::bc_elem1
 * @param alpha
 * @param ell
 * @param beta
 * @return
 */
int SemiHomopolymerAlignmentStat::bc_elem1(int strand, char alpha, int ell, char beta){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE)) |
                 (ell<<HomopolymerSpace::BETABITSIZE) | (HomopolymerSpace::index(beta));
    return base_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::qc_elem1
 * @param alpha
 * @param ell
 * @param qual
 * @return
 */
int SemiHomopolymerAlignmentStat::qc_elem1(int strand, char alpha, int ell, int qual){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE)) |
                 (ell<<HomopolymerSpace::QUALITYBITSIZE) | qual;
    return qual_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::lc_elem1
 * @param alpha
 * @param ell
 * @param kappa
 * @return
 */
int SemiHomopolymerAlignmentStat::lc_elem1(int strand, char alpha, int ell, int kappa){
    uint32_t idx=(strand<<(HomopolymerSpace::ALPHABITSIZE+HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE) |
                  HomopolymerSpace::index(alpha)<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE)) |
                 (ell<<HomopolymerSpace::LENGTHBITSIZE) | kappa;
    return len_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::print
 */
void SemiHomopolymerAlignmentStat::print(){
    char buffer[2048];
    char alpha,beta;

    for (int nucl=HomopolymerSpace::A; nucl<HomopolymerSpace::NucleotideSize; nucl++){
        alpha=HomopolymerSpace::symbol((HomopolymerSpace::Nucleotide)nucl)[0];
        for (int strand=0; strand<=1; strand++){
            string direction=(strand==1)?"Positive":"Negative";
            // general title
            cout<<"[Statistics of Homopolymer-Nucleotide "<<alpha<<" Strand:"<<direction<<"]"<<endl;
            // base-call table
            cout<<"[#1: Base-Call Table]"<<endl;
            sprintf(buffer,"      %10s %10s %10s %10s","A","C","G","T");
            cout<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int nucl2=HomopolymerSpace::A; nucl2<=HomopolymerSpace::T; nucl2++){
                    beta=HomopolymerSpace::symbol((HomopolymerSpace::Nucleotide)nucl2)[0];
                    sprintf(buffer,"%s %10d",buffer,bc_elem1(strand,alpha,ell,beta));
                }
                cout<<buffer<<endl;
            }
            // qual-call table
            cout<<"[#2: Quality-Score-Call Table]"<<endl;
            sprintf(buffer,"size ");
            for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual++){
                sprintf(buffer, "%s %10d", buffer, qual);
            }
            cout<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual++){
                    sprintf(buffer, "%s %10d",buffer,qc_elem1(strand,alpha,ell,qual));
                }
                cout<<buffer<<endl;
            }
            // length-call table
            cout<<"[#3: Length-Call Table]"<<endl;
            sprintf(buffer,"size ");
            for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                sprintf(buffer,"%s %10d",buffer, kappa);
            }
            cout<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                    sprintf(buffer, "%s %10d",buffer,lc_elem1(strand,alpha,ell,kappa));
                }
                cout<<buffer<<endl;
            }
            cout<<endl<<endl;
        }
    }
}

/**
 * @brief SemiHomopolymerAlignmentStat::print
 * @param filename
 */
void SemiHomopolymerAlignmentStat::print(string filename){
    char buffer[2048];
    char alpha,beta;
    for (int nucl=HomopolymerSpace::A; nucl<HomopolymerSpace::NucleotideSize; nucl++){
        alpha=HomopolymerSpace::symbol((HomopolymerSpace::Nucleotide)nucl)[0];
        for (int strand=0; strand<=1; strand++){
            string direction=(strand==1)?"positive":"negative";

            string fn=filename+"."+HomopolymerSpace::symbol((HomopolymerSpace::Nucleotide)nucl)+"."+direction;
            ofstream ofs(fn.c_str());

            // base-call table
            ofs<<"[#1: Base-Call Table]"<<endl;
            sprintf(buffer,"Size  %10s %10s %10s %10s","A","C","G","T");
            ofs<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int nucl2=HomopolymerSpace::A; nucl2<=HomopolymerSpace::T; nucl2++){
                    beta=HomopolymerSpace::symbol((HomopolymerSpace::Nucleotide)nucl2)[0];
                    sprintf(buffer,"%s %10d",buffer,bc_elem1(strand,alpha,ell,beta));
                }
                ofs<<buffer<<endl;
            }
            // qual-call table
            ofs<<"[#2: Quality-Score-Call Table]"<<endl;
            sprintf(buffer,"size ");
            for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual++){
                sprintf(buffer, "%s %10d", buffer, qual);
            }
            ofs<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual++){
                    sprintf(buffer, "%s %10d",buffer,qc_elem1(strand,alpha,ell,qual));
                }
                ofs<<buffer<<endl;
            }
            // length-call table
            ofs<<"[#3: Length-Call Table]"<<endl;
            sprintf(buffer,"size ");
            for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                sprintf(buffer,"%s %10d",buffer, kappa);
            }
            ofs<<buffer<<endl;
            for (int ell=0; ell<HomopolymerSpace::HomopolymerSizeMax; ell++){
                sprintf(buffer,"%4d ",ell);
                for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                    sprintf(buffer, "%s %10d",buffer,lc_elem1(strand,alpha,ell,kappa));
                }
                ofs<<buffer<<endl;
            }
            ofs.close();
        }
    }
}



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
    align_query=homo_align.align_query;
    align_target=homo_align.align_target;
    align_status=homo_align.align_status;

    this->align_name=homo_align.align_name;
    this->query_seq=homo_align.query_seq;
    this->query_qual=homo_align.query_qual;
    this->query_strand=homo_align.query_strand;
    this->target_seq=homo_align.target_seq;
}

/**
 * @brief SemiHomopolymerAlignment::SemiHomopolymerAlignment
 * @param nucl_align
 */
SemiHomopolymerAlignment::SemiHomopolymerAlignment(const NucleotideAlignment &nucl_align){
    this->align_name=nucl_align.align_name;
    this->query_seq=nucl_align.query_seq;
    this->query_qual=nucl_align.query_qual;
    this->query_strand=nucl_align.query_strand;
    this->target_seq=nucl_align.target_seq;
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
    align_query=homo_align.align_query;
    align_target=homo_align.align_target;
    align_status=homo_align.align_status;

    this->align_name=homo_align.align_name;
    this->query_seq=homo_align.query_seq;
    this->query_qual=homo_align.query_qual;
    this->query_strand=homo_align.query_strand;
    this->target_seq=homo_align.target_seq;

    return *this;
}

SemiHomopolymerAlignment& SemiHomopolymerAlignment::operator =(const NucleotideAlignment& nucl_align){
    this->align_name=nucl_align.align_name;
    this->query_seq=nucl_align.query_seq;
    this->query_qual=nucl_align.query_qual;
    this->query_strand=nucl_align.query_strand;
    this->target_seq=nucl_align.target_seq;

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
 * @brief SemiHomopolymerAlignment::statistics
 * @param stat
 * @param cycles
 */
void SemiHomopolymerAlignment::statistics(vector<SemiHomopolymerAlignmentStat> &stat, int cycles){
    // TODO: consider the alignment status
    // TODO: consider higher order markov statistics

    // clear everything
    stat=vector<SemiHomopolymerAlignmentStat>(cycles);

    // target length without spaces
    int len=0;
    for (int i=0; i<align_target.len; i++){
        len+=align_target.ell[i];
    }

    // scan through the alignment
    for (int i=0,j=0; i<align_target.len; i++,j+=align_target.ell[i]){
        // cycle no.
        int c=(int)(j*cycles/(len+0.));
        // statistician
        SemiHomopolymerAlignmentStat s=stat[c];
        // symbols and size
        char alpha,beta;
        int ell,kappa,qual;

        alpha=align_target.alpha[i][0];
        ell=align_target.ell[i];
        //if target is a space at current position
        if (align_target.ell[i]==0){
            alpha='-';
        }

        if (HomopolymerSpace::index(alpha)<HomopolymerSpace::NucleotideSize){
            kappa=align_query.seq[i].size();
            for (int k=0; k<kappa; k++){
                beta=align_query.seq[i][k];
                qual=(int)align_query.qual[i][k]-33;
                if (HomopolymerSpace::index(beta)<HomopolymerSpace::NucleotideSize){
                    // update base-call table
                    s.bc_incr1(query_strand, alpha, ell, beta);
                    // update quality-score-call table
                    s.qc_incr1(query_strand, alpha, ell, qual);
                }
            }
            // update length-call table
            s.lc_incr1(query_strand, alpha, ell, kappa);
        }

        // resotre the statistics
        stat[c]=s;
    }
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

/**
 * @brief SemiHomopolymerAlignmentPool::statistics
 * @param cycles
 */
void SemiHomopolymerAlignmentPool::statistics(int cycles){
   vector<SemiHomopolymerAlignmentStat> stat;
   // do the statistics
   statistics(stat,cycles);
   // print out the statistics
   for (int c=0; c<cycles; c++){
       cout<<"[Cycle #"<<c<<"]"<<endl;
       stat[c].print();
   }
}

void SemiHomopolymerAlignmentPool::statistics(string filename, int cycles){
    vector<SemiHomopolymerAlignmentStat> stat;
    // do the statistics
    statistics(stat,cycles);
    // print out the statistics
    for (int c=0; c<cycles; c++){
        stat[c].print(filename+".c"+to_string(c));
    }

}

/**
 * @brief SemiHomopolymerAlignmentPool::statistics
 * @param stat
 * @param cycles
 */
void SemiHomopolymerAlignmentPool::statistics(vector<SemiHomopolymerAlignmentStat> &stat, int cycles){
    // clear everything
    stat=vector<SemiHomopolymerAlignmentStat>(cycles);

    // loop over the alignments
    for (int i=0; i<(int)align_pool.size(); i++){
        // current alignment
        SemiHomopolymerAlignment aln1=align_pool[i];
        // statistician
        vector<SemiHomopolymerAlignmentStat> aln1_stat;
        // do the statistics
        aln1.statistics(aln1_stat,cycles);
        // update the overall statistician
        for (int c=0; c<cycles; c++){
            stat[c]+=aln1_stat[c];
        }
    }
}
