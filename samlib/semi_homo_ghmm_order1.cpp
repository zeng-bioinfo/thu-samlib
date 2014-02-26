#include "semi_homo_ghmm_order1.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "glm.h"
#include <iostream>
#include <fstream>
using namespace std;

/**
 * @brief SemiHomopolymerGHMMOrder1::SemiHomopolymerGHMMOrder1
 */
SemiHomopolymerGHMMOrder1::SemiHomopolymerGHMMOrder1(){
    // default the read is a segment
    this->cycles=1;
    initialization();
}

/**
 * @brief SemiHomopolymerGHMMOrder1::SemiHomopolymerGHMMOrder1
 * @param _cycles
 */
SemiHomopolymerGHMMOrder1::SemiHomopolymerGHMMOrder1(int _cycles){
    this->cycles=_cycles;
    initialization();
}

/**
 * @brief SemiHomopolymerGHMMOrder1::~SemiHomopolymerGHMMOrder1
 */
SemiHomopolymerGHMMOrder1::~SemiHomopolymerGHMMOrder1(){
    delete [] this->log_prob_state_trans;
    delete [] this->log_prob_base_call;
    delete [] this->log_prob_qual_call;
    delete [] this->log_prob_length_call;
    delete [] this->glm_poisson_b0;
    delete [] this->glm_poisson_b1;
}

/**
 * @brief SemiHomopolymerGHMMOrder1::initialization
 */
void SemiHomopolymerGHMMOrder1::initialization(){
    size_prob_state_trans_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
            SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE;
    size_prob_base_call_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
            HomopolymerSpace::HomopolymerPosMax*
            HomopolymerSpace::QualityScoreSlot*
            NucleotideSpace::NucleotideSize;
    size_prob_qual_call_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
            HomopolymerSpace::HomopolymerPosMax*
            NucleotideSpace::QualityScoreSlot;
    size_prob_length_call_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
            SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX;
    size_glm_poisson_b0_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            NucleotideSpace::NucleotideSize;
    size_glm_poisson_b1_table=(1<<HomopolymerSpace::STRANDBITSIZE)*cycles*
            NucleotideSpace::NucleotideSize;

    this->log_prob_state_trans=new double[size_prob_state_trans_table];
    this->log_prob_base_call=new double[size_prob_base_call_table];
    this->log_prob_length_call=new double[size_prob_length_call_table];
    this->log_prob_qual_call=new double[size_prob_qual_call_table];
    this->glm_poisson_b0=new double[size_glm_poisson_b0_table];
    this->glm_poisson_b1=new double[size_glm_poisson_b1_table];

    // set to zero
    memset(this->log_prob_state_trans, 0, sizeof(double)*size_prob_state_trans_table);
    memset(this->log_prob_base_call, 0, sizeof(double)*size_prob_base_call_table);
    memset(this->log_prob_length_call, 0, sizeof(double)*size_prob_length_call_table);
    memset(this->log_prob_qual_call, 0, sizeof(double)*size_prob_qual_call_table);
    memset(this->glm_poisson_b0, 0, sizeof(double)*size_glm_poisson_b0_table);
    memset(this->glm_poisson_b1, 0, sizeof(double)*size_glm_poisson_b1_table);
}

/**
 * @brief SemiHomopolymerGHMMOrder1::reset
 */
void SemiHomopolymerGHMMOrder1::reset(){
    // set to zero
    memset(this->log_prob_state_trans, 0, sizeof(double)*size_prob_state_trans_table);
    memset(this->log_prob_base_call, 0, sizeof(double)*size_prob_base_call_table);
    memset(this->log_prob_length_call, 0, sizeof(double)*size_prob_length_call_table);
    memset(this->log_prob_qual_call, 0, sizeof(double)*size_prob_qual_call_table);
    memset(this->glm_poisson_b0, 0, sizeof(double)*size_glm_poisson_b0_table);
    memset(this->glm_poisson_b1, 0, sizeof(double)*size_glm_poisson_b1_table);
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_estimate
 * @param align_pool
 * @param maxIter
 * @param thres
 */
void SemiHomopolymerGHMMOrder1::parameter_estimate(SemiHomopolymerAlignmentPool &align_pool, int band, int maxIter, double thres){
    cout << "Train the parameters of HMM" << endl;

    int iter;
    double lls0,lls1;   // lls0 is sum of log-likelihood-scores at previous time
                        // lls1 is sum of log-likelihood-scores at current time
    double dlls;
    // items for temporary use
    double *tmp_log_prob_base_call; // temporary storage of log_prob_base_call
    double *tmp_log_prob_qual_call; // temporary storage of log_prob_qual_call
    double *tmp_log_prob_state_trans;   // temporary storage of hidden state transition
    double *tmp_glm_poisson_b0;     // temporary storage of poisson parameter
    double *tmp_glm_poisson_b1;     // temporary storage of poisson parameter

    tmp_log_prob_state_trans=new double[size_prob_state_trans_table];
    tmp_log_prob_base_call=new double[size_prob_base_call_table];
    tmp_log_prob_qual_call=new double[size_prob_qual_call_table];
    tmp_glm_poisson_b0=new double[size_glm_poisson_b0_table];
    tmp_glm_poisson_b1=new double[size_glm_poisson_b1_table];

    // TODO: compute the initial parameters from given data
    // debug
    parameter_initialize(align_pool);

    // TODO: iteratively update the parameters
    // while (iter<maxIter && lls1-lls0>thres)
    //      save old parameters into temporary storage
    //      compute new alignments using old parameters
    //      compute log-likelihood-score of old parameters
    //      compute new parameters from new alignments
    //      compute log-likelihood-score of new parameters

    vector<double> every_align_change(align_pool.align_pool.size(),100);
    vector<double> score0(align_pool.align_pool.size(),-1e+300);
    vector<SemiHomopolymerAlignment> align0(align_pool.align_pool.size());

    vector<SemiHomopolymerAlignmentStat> statistics(cycles);


    lls0=0;
    for (int i=0; i<(int)align_pool.align_pool.size(); i++){
        lls0+=compute_log_likelihood_score(align_pool.align_pool[i]);
    }
    if (lls0!=lls0){
        cout<<"Bad initial parameters\n"
            <<"Please run again\n"
            <<"Exit"<<endl;
        exit(0);
    }

    iter=0;
    lls1=1e+300;
    dlls=1e+300;
    while (iter<maxIter && fabs(dlls)>thres) {
        cout << "[Iteration #" << iter << "]" << endl;
        // temporary storage
        memcpy(tmp_log_prob_state_trans, this->log_prob_state_trans, size_prob_state_trans_table*sizeof(double));
        memcpy(tmp_log_prob_base_call, this->log_prob_base_call, size_prob_base_call_table*sizeof(double));
        memcpy(tmp_log_prob_qual_call, this->log_prob_qual_call, size_prob_qual_call_table*sizeof(double));
        memcpy(tmp_glm_poisson_b0, this->glm_poisson_b0, size_glm_poisson_b0_table*sizeof(double));
        memcpy(tmp_glm_poisson_b1, this->glm_poisson_b1, size_glm_poisson_b1_table*sizeof(double));

        lls1=0;
        // compute new alignments using old parameters
        for (int i=0; i<(int)align_pool.align_pool.size(); i++){
            {
                char buffer[1024];
                // show progress
                int percent=(int)((i+1)*100.0/(align_pool.align_pool.size()+0.));

                string bar;

                for(int ii = 0; ii < 50; ii++){
                    if( ii < (percent/2)){
                        bar.replace(ii,1,"=");
                    }else if( ii == (percent/2)){
                        bar.replace(ii,1,">");
                    }else{
                        bar.replace(ii,1," ");
                    }
                }

                cout<< "\r";
                cout<< "Compute " << align_pool.align_pool.size() << " re-alignments   \t"
                    << "[" << bar << "] ";
                sprintf(buffer, "%.2f", ((i+1)*100.0/(align_pool.align_pool.size()+0.)));
                cout<< buffer << "%     ";

                cout << flush;
            }

            if (every_align_change[i]>0.01){
                SemiHomopolymerAlignment tmp_aln;
                double tmp_score=compute_banded_realignment(align_pool.align_pool[i], tmp_aln, band);
                if (tmp_score>score0[i]){
                    lls1+=tmp_score;
                    align_pool.align_pool[i]=tmp_aln;

                    every_align_change[i]=tmp_score-score0[i];
                    score0[i]=tmp_score;
                    align0[i]=tmp_aln;
                }else{
                    lls1+=score0[i];
                    align_pool.align_pool[i]=align0[i];
                    every_align_change[i]=0;
                }
            }else{
                lls1+=score0[i];
                align_pool.align_pool[i]=align0[i];
            }
        }

        {
            // show progress
            cout<<endl;
        }

        // counting
        align_pool.statistics(statistics, cycles);

        // compute new parameters
        cout<<"Compute new parameters"<<endl;
        parameter_update(statistics);

        // change of log-likelihood score
        dlls=lls1-lls0;

        // iteration log
        cout.precision(15);
        cout<<"lls("<<iter-1<<")="<<lls0
            <<" lls("<<iter<<")="<<lls1
            <<" lls("<<iter<<")-lls("<<iter-1<<")="<<dlls<<endl;

        lls0=lls1;

        iter++;
    }

    if (dlls<0){
        memcpy(this->log_prob_state_trans, tmp_log_prob_state_trans, size_prob_state_trans_table*sizeof(double));
        memcpy(this->log_prob_base_call, tmp_log_prob_base_call, size_prob_base_call_table*sizeof(double));
        memcpy(this->log_prob_qual_call, tmp_log_prob_qual_call, size_prob_qual_call_table*sizeof(double));
        memcpy(this->glm_poisson_b0, tmp_glm_poisson_b0, size_glm_poisson_b0_table*sizeof(double));
        memcpy(this->glm_poisson_b1, tmp_glm_poisson_b1, size_glm_poisson_b1_table*sizeof(double));
    }

    cout<<"Complete the model training"<<endl;

    delete [] tmp_log_prob_base_call;
    delete [] tmp_log_prob_qual_call;
    delete [] tmp_log_prob_state_trans;
    delete [] tmp_glm_poisson_b0;
    delete [] tmp_glm_poisson_b1;
    vector<double>().swap(every_align_change);
    vector<double>().swap(score0);
    vector<SemiHomopolymerAlignment>().swap(align0);
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_initialize
 * @param align_pool
 */
void SemiHomopolymerGHMMOrder1::parameter_initialize(SemiHomopolymerAlignmentPool &align_pool){
    // #1: count the occurrences of patterns from the given alignments
    vector<SemiHomopolymerAlignmentStat> stat;
    align_pool.statistics(stat, cycles);

    // #2: update the parameters
    cout << "Compute initial parameters" << endl;
    parameter_update(stat);
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_update
 * @param statistics
 */
void SemiHomopolymerGHMMOrder1::parameter_update(vector<SemiHomopolymerAlignmentStat> &statistics){
    // clear old parameters
    reset();

    // compute glm parameters using overall data
    int idx0,idx1;
    int ell;

    vector<double> pseudo(20*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX,0);
    for (int c=0; c<cycles; c++){
        for (int strand=0; strand<=1; strand++){
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
                ell=SemiHomopolymerAlignmentSpace::state_length[pi];
                for (int k=0; k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                         pi*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                         k;
                    idx1=(ell-1)*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                         k;
                    pseudo[idx1]+=statistics[c].len_call_table[idx0];
                }
            }
        }
    }

    vector<double> X, y;
    double gb[2];
    int n=0;
    for (ell=1; ell<=20; ell++){
        idx1=(ell-1)*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
             ell;
        if (pseudo[idx1]>100){
            int dk;
            double tot=0,Z=0;
            for (int k=0; k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                dk=abs(ell-k);
                idx1=(ell-1)*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                     k;
                tot+=dk*pseudo[idx1];
                Z+=pseudo[idx1];
            }
            X.push_back(ell);
            y.push_back(tot/Z);
            n++;
        }
    }
    GLM::irls(X.data(),y.data(),n,gb,100,1e-7);

    if (gb[0]!=gb[0] || gb[1]!=gb[1]){
        gb[0]=log(1e-5);
        gb[1]=0;
    }


    // update parameters
    for (int c=0; c<cycles; c++){
        // #1: state transition probability
        parameter_update_state_trans_prob(statistics[c].state_trans_table, c);
        // #2: base calling probability
        parameter_update_base_call_prob(statistics[c].base_call_table, c);
        // #3: quality calling probability
        parameter_update_qual_call_prob(statistics[c].qual_call_table, c);
        // #4: length calling probability
        parameter_update_len_call_prob(statistics[c].len_call_table, c, gb);
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_update_state_trans_prob
 * @param st_count
 * @param c
 */
void SemiHomopolymerGHMMOrder1::parameter_update_state_trans_prob(const vector<int> &st_count, int c){
    int idx;
    int idx0;

    // plus the count and pseudo-count
    for (int strand=1; strand>=0; strand--){
        // 1. B
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ((int)SemiHomopolymerAlignmentSpace::B)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                (int)SemiHomopolymerAlignmentSpace::SH;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::B);
        log_prob_state_trans[idx]=st_count[idx0]+.1;
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ((int)SemiHomopolymerAlignmentSpace::B)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                pi;
            idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::B);
            if (pi==(int)SemiHomopolymerAlignmentSpace::A1 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::C1 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::G1 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::T1){
                log_prob_state_trans[idx]=st_count[idx0]+100;
            }else{
                log_prob_state_trans[idx]=st_count[idx0]+.1;
            }
        }

        // 2. SH
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            ((int)SemiHomopolymerAlignmentSpace::SH)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            (int)SemiHomopolymerAlignmentSpace::SH;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::SH);
        log_prob_state_trans[idx]=st_count[idx0]+100;
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ((int)SemiHomopolymerAlignmentSpace::SH)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                pi;
            idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::SH);
            log_prob_state_trans[idx]=st_count[idx0]+.1;
        }

        // 3. from A1..T20
        for (int ppi=(int)SemiHomopolymerAlignmentSpace::A1; ppi<=(int)SemiHomopolymerAlignmentSpace::T20; ppi++){
            // to A1..T20
            for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
                if ((ppi>=(int)SemiHomopolymerAlignmentSpace::A1 && ppi<=(int)SemiHomopolymerAlignmentSpace::A20) &&
                        (pi>=(int)SemiHomopolymerAlignmentSpace::A1 && pi<=(int)SemiHomopolymerAlignmentSpace::A20)) continue;
                if ((ppi>=(int)SemiHomopolymerAlignmentSpace::C1 && ppi<=(int)SemiHomopolymerAlignmentSpace::C20) &&
                        (pi>=(int)SemiHomopolymerAlignmentSpace::C1 && pi<=(int)SemiHomopolymerAlignmentSpace::C20)) continue;
                if ((ppi>=(int)SemiHomopolymerAlignmentSpace::G1 && ppi<=(int)SemiHomopolymerAlignmentSpace::G20) &&
                        (pi>=(int)SemiHomopolymerAlignmentSpace::G1 && pi<=(int)SemiHomopolymerAlignmentSpace::G20)) continue;
                if ((ppi>=(int)SemiHomopolymerAlignmentSpace::T1 && ppi<=(int)SemiHomopolymerAlignmentSpace::T20) &&
                        (pi>=(int)SemiHomopolymerAlignmentSpace::T1 && pi<=(int)SemiHomopolymerAlignmentSpace::T20)) continue;
                idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                    ppi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                    pi;
                idx=state_trans_table_index(strand, c, pi, ppi);
                log_prob_state_trans[idx]=st_count[idx0]+.1;
            }
            // to I
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ppi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                (int)SemiHomopolymerAlignmentSpace::I;
            idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::I, ppi);
            log_prob_state_trans[idx]=st_count[idx0]+.1;
            // to ST
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ppi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                (int)SemiHomopolymerAlignmentSpace::ST;
            idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::ST, ppi);
            log_prob_state_trans[idx]=st_count[idx0]+.1;
            // to E
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ppi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                (int)SemiHomopolymerAlignmentSpace::E;
            idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::E, ppi);
            log_prob_state_trans[idx]=st_count[idx0]+.1;
        }

        // 4. I
        // to A1..T20
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ((int)SemiHomopolymerAlignmentSpace::I)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                pi;
            idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::I);
            log_prob_state_trans[idx]=st_count[idx0]+.1;
        }
        // to I
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            ((int)SemiHomopolymerAlignmentSpace::I)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            (int)SemiHomopolymerAlignmentSpace::I;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::I, (int)SemiHomopolymerAlignmentSpace::I);
        log_prob_state_trans[idx]=st_count[idx0]+.1;
        // to ST
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            ((int)SemiHomopolymerAlignmentSpace::I)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            (int)SemiHomopolymerAlignmentSpace::ST;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::ST, (int)SemiHomopolymerAlignmentSpace::I);
        log_prob_state_trans[idx]=st_count[idx0]+.1;

        // 5. ST
        // to ST
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            ((int)SemiHomopolymerAlignmentSpace::ST)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            (int)SemiHomopolymerAlignmentSpace::ST;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::ST, (int)SemiHomopolymerAlignmentSpace::ST);
        log_prob_state_trans[idx]=st_count[idx0]+1;
        // to E
        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            ((int)SemiHomopolymerAlignmentSpace::ST)*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
            (int)SemiHomopolymerAlignmentSpace::E;
        idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::E, (int)SemiHomopolymerAlignmentSpace::ST);
        log_prob_state_trans[idx]=st_count[idx0]+.1;

        // compute the probability
        for (int ppi=(int)SemiHomopolymerAlignmentSpace::A1; ppi<(int)SemiHomopolymerAlignmentSpace::E; ppi++){
            double norm=0;
            for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                idx=state_trans_table_index(strand, c, pi, ppi);
                norm+=log_prob_state_trans[idx];
            }
            for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                idx=state_trans_table_index(strand, c, pi, ppi);
                if (log_prob_state_trans[idx]>0){
                    log_prob_state_trans[idx]=log(log_prob_state_trans[idx])-log(norm+1e-7);
                }else{
                    log_prob_state_trans[idx]=LOGZERO;
                }
            }
        }
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_update_base_call_prob
 * @param bc_count
 * @param c
 */
void SemiHomopolymerGHMMOrder1::parameter_update_base_call_prob(const vector<int> &bc_count, int c){
    int idx,idx0;
    int last_idx00, last_idx01;
    int last_idx10, last_idx11;
    int last_idx20, last_idx21;
    vector<int> count(NucleotideSpace::NucleotideSize, 0);
    int num,max_num;
    double norm;

    for (int strand=1; strand>=0; strand--){
        for (int PI=(int)SemiHomopolymerAlignmentSpace::A1; PI<=(int)SemiHomopolymerAlignmentSpace::T20; PI+=20){
            // compute pseudo count
            count.assign(count.size(), 0);
            for (int pi=PI; pi<PI+20; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    q*NucleotideSpace::NucleotideSize+
                                    beta;
                            count[beta]+=bc_count[idx0];
                        }
                    }
                }
            }
            for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                if (NucleotideSpace::idx2nucl[beta]==SemiHomopolymerAlignmentSpace::state_alpha[PI]){
                    count[beta]+=1000;
                }else{
                    count[beta]+=1;
                }
            }

            // level state
            last_idx00=-1;
            last_idx01=-2;
            for (int pi=PI; pi<PI+20; pi++){
                max_num=0;
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        num=0;
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    q*NucleotideSpace::NucleotideSize+
                                    beta;
                            num+=bc_count[idx0];
                        }
                        if (max_num<num) max_num=num;
                    }
                }
                if (last_idx00<0 && max_num>=count_thresh) last_idx00=pi;
                if (max_num>=count_thresh) last_idx01=pi;
            }
            // level position
            // [last_idx00, last_idx01]
            for (int pi=last_idx00; pi<=last_idx01; pi++){
                last_idx10=-1;
                last_idx11=-2;
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    max_num=0;
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        num=0;
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    q*NucleotideSpace::NucleotideSize+
                                    beta;
                            num+=bc_count[idx0];
                        }
                        if (max_num<num) max_num=num;
                    }
                    if (last_idx10<0 && max_num>=count_thresh) last_idx10=i;
                    if (max_num>=count_thresh) last_idx11=i;
                }
                // [last_idx10, last_idx11]
                for (int i=last_idx10; i<=last_idx11; i++){
                    // level quality
                    last_idx20=-1;
                    last_idx21=-2;
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        num=0;
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    q*NucleotideSpace::NucleotideSize+
                                    beta;
                            num+=bc_count[idx0];
                        }
                        if (last_idx20<0 && num>=count_thresh) last_idx20=q;
                        if (num>=count_thresh) last_idx21=q;
                    }
                    // [last_idx20,last_idx21]
                    for (int q=last_idx20; q<=last_idx21; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                                    q*NucleotideSpace::NucleotideSize+
                                    beta;
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            this->log_prob_base_call[idx]=bc_count[idx0]+1;
                        }
                    }
                    // [0,last_idx20)
                    for (int q=0; q<last_idx20; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=this->base_call_table_index(strand, c, pi, i, last_idx20, beta);
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                        }
                    }
                    // (last_idx21,inf]
                    last_idx21=(last_idx21<0)?-1:last_idx21;
                    for (int q=last_idx21+1; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            if (last_idx21==-1){
                                this->log_prob_base_call[idx]=count[beta];
                            }else{
                                idx0=this->base_call_table_index(strand, c, pi, i, last_idx21, beta);
                                this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                            }
                        }
                    }
                }
                // [0,last_idx10)
                for (int i=0; i<last_idx10; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=this->base_call_table_index(strand, c, pi, last_idx10, q, beta);
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                        }
                    }
                }
                // (last_idx11,inf]
                last_idx11=(last_idx11<0)?-1:last_idx11;
                for (int i=last_idx11+1; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            if (last_idx11==-1){
                                this->log_prob_base_call[idx]=count[beta];
                            }else{
                                idx0=this->base_call_table_index(strand, c, pi, last_idx11, q, beta);
                                this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                            }
                        }
                    }
                }
            }
            // [0,last_00)
            for (int pi=PI; pi<last_idx00; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx0=this->base_call_table_index(strand, c, last_idx00, i, q, beta);
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                        }
                    }
                }
            }
            // (last_01,inf]
            last_idx01=(last_idx01<0)?PI-1:last_idx01;
            for (int pi=last_idx01+1; pi<PI+20; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                            if (last_idx01==PI-1){
                                this->log_prob_base_call[idx]=count[beta];
                            }else{
                                idx0=this->base_call_table_index(strand, c, last_idx01, i, q, beta);
                                this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                            }
                        }
                    }
                }
            }
        }

        // <state: I, SH, ST>
        // compute last_idx0 for level pi
        for (int pi=(int)SemiHomopolymerAlignmentSpace::I; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
            // compute pseudo count
            count.assign(count.size(), 0);
            int i=0;
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            q*NucleotideSpace::NucleotideSize+
                            beta;
                    count[beta]+=bc_count[idx0];
                }
            }
            for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                count[beta]+=1;
            }
            // compute table
            // level quality
            last_idx20=-1;
            last_idx21=-2;
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                num=0;
                for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            q*NucleotideSpace::NucleotideSize+
                            beta;
                    num+=bc_count[idx0];
                }
                if (last_idx20<0 && num>=count_thresh) last_idx20=q;
                if (num>=count_thresh) last_idx21=q;
            }
            // [last_idx20,last_idx21]
            for (int q=last_idx20; q<=last_idx21; q++){
                for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            i*NucleotideSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                            q*NucleotideSpace::NucleotideSize+
                            beta;
                    idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                    this->log_prob_base_call[idx]=bc_count[idx0]+1;
                }
            }
            // [0,last_idx20)
            for (int q=0; q<last_idx20; q++){
                for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                    idx0=this->base_call_table_index(strand, c, pi, i, last_idx20, beta);
                    idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                    this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                }
            }
            // (last_idx21,inf]
            last_idx21=(last_idx21<0)?-1:last_idx21;
            for (int q=last_idx21+1; q<NucleotideSpace::QualityScoreSlot; q++){
                for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                    idx=this->base_call_table_index(strand, c, pi, i, q, beta);
                    if (last_idx21==-1){
                        this->log_prob_base_call[idx]=count[beta];
                    }else{
                        idx0=this->base_call_table_index(strand, c, pi, i, last_idx21, beta);
                        this->log_prob_base_call[idx]=this->log_prob_base_call[idx0];
                    }
                }
            }
        }

        // compute the probability
        for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
            for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                    norm=0;
                    for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                        idx=base_call_table_index(strand, c, pi, i, q, beta);
                        norm+=log_prob_base_call[idx];
                    }
                    for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                        idx=base_call_table_index(strand, c, pi, i, q, beta);
                        if (log_prob_base_call[idx]>0){
                            log_prob_base_call[idx]=log(log_prob_base_call[idx])-log(norm+1e-7);
                        }else{
                            log_prob_base_call[idx]=LOGZERO;
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_update_qual_call_prob
 * @param qc_count
 * @param c
 */
void SemiHomopolymerGHMMOrder1::parameter_update_qual_call_prob(const vector<int> &qc_count, int c){
    for (int strand=1; strand>=0; strand--){
        int idx,idx0;
        int last_idx00, last_idx01;
        int last_idx10, last_idx11;
        int num,max_num;
        double norm;
        vector<int> count(NucleotideSpace::QualityScoreSlot);
        // compute the table
        // state A1..T20
        for (int PI=(int)SemiHomopolymerAlignmentSpace::A1; PI<=(int)SemiHomopolymerAlignmentSpace::T20; PI+=20){
            // compute pseudo count
            count.assign(count.size(), 0);
            for (int pi=PI; pi<PI+20; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                i*NucleotideSpace::QualityScoreSlot+
                                q;
                        count[q]+=qc_count[idx0];
                    }
                }
            }
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                count[q]+=1;
            }
            // level state
            last_idx00=-1;
            last_idx01=-2;
            for (int pi=PI; pi<PI+20; pi++){
                max_num=0;
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    num=0;
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                i*NucleotideSpace::QualityScoreSlot+
                                q;
                        num+=qc_count[idx0];
                    }
                    if (max_num<num) max_num=num;
                }
                if (last_idx00<0 && max_num>=count_thresh) last_idx00=pi;
                if (max_num>=count_thresh) last_idx01=pi;
            }
            // [last_idx00,last_idx01]
            for (int pi=last_idx00; pi<=last_idx01; pi++){
                // level position
                last_idx10=-1;
                last_idx11=-2;
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    num=0;
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                i*NucleotideSpace::QualityScoreSlot+
                                q;
                        num+=qc_count[idx0];
                    }
                    if (last_idx10<0 && num>=count_thresh) last_idx10=i;
                    if (num>=count_thresh) last_idx11=i;
                }
                // [last_idx10,last_idx11]
                for (int i=last_idx10; i<=last_idx11; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                                i*NucleotideSpace::QualityScoreSlot+
                                q;
                        idx=this->quality_call_table_index(strand, c, pi, i, q);
                        this->log_prob_qual_call[idx]=qc_count[idx0]+1;
                    }
                }
                // [0,last_idx10)
                for (int i=0; i<last_idx10; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=this->quality_call_table_index(strand, c, pi, last_idx10, q);
                        idx=this->quality_call_table_index(strand, c, pi, i, q);
                        this->log_prob_qual_call[idx]=this->log_prob_qual_call[idx0];
                    }
                }
                // (last_idx11,inf]
                last_idx11=(last_idx11<0)?-1:last_idx11;
                for (int i=last_idx11+1; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx=this->quality_call_table_index(strand, c, pi, i, q);
                        if (last_idx11==-1){
                            this->log_prob_qual_call[idx]=count[q];
                        }else{
                            idx0=this->quality_call_table_index(strand, c, pi, last_idx11, q);
                            this->log_prob_qual_call[idx]=this->log_prob_qual_call[idx0];
                        }
                    }
                }
            }
            // [0,last_idx00)
            for (int pi=PI; pi<last_idx00; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx0=this->quality_call_table_index(strand, c, last_idx00, i, q);
                        idx=this->quality_call_table_index(strand, c, pi, i, q);
                        this->log_prob_qual_call[idx]=this->log_prob_qual_call[idx0];
                    }
                }
            }
            // (last_idx01,inf]
            last_idx01=(last_idx01<0)?PI-1:last_idx01;
            for (int pi=last_idx01+1; pi<PI+20; pi++){
                for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                    for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                        idx=this->quality_call_table_index(strand, c, pi, i, q);
                        if (last_idx01==PI-1){
                            this->log_prob_qual_call[idx]=count[q];
                        }else{
                            idx0=this->quality_call_table_index(strand, c, last_idx01, i, q);
                            this->log_prob_qual_call[idx]=this->log_prob_qual_call[idx0];
                        }
                    }
                }
            }
        }
        // state I,SH,ST
        for (int pi=(int)SemiHomopolymerAlignmentSpace::I; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
            // compute pseudo count
            count.assign(count.size(), 0);
            int i=0;
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                        pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                        i*NucleotideSpace::QualityScoreSlot+
                        q;
                count[q]+=qc_count[idx0];
            }
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                count[q]+=1;
            }
            // compute table
            for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                idx=this->quality_call_table_index(strand, c, pi, i, q);
                this->log_prob_qual_call[idx]=count[q];
            }
        }

        // compute the probability
        for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
            for (int i=0; i<HomopolymerSpace::HomopolymerPosMax; i++){
                norm=0;
                for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                    idx=quality_call_table_index(strand, c, pi, i, q);
                    norm+=log_prob_qual_call[idx];
                }
                for (int q=0; q<NucleotideSpace::QualityScoreSlot; q++){
                    idx=quality_call_table_index(strand, c, pi, i, q);
                    if (log_prob_qual_call[idx]>0){
                        log_prob_qual_call[idx]=log(log_prob_qual_call[idx])-log(norm+1e-7);
                    }else{
                        log_prob_qual_call[idx]=LOGZERO;
                    }
                }
            }
        }
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_update_len_call_prob
 * @param lc_count
 * @param c
 */
void SemiHomopolymerGHMMOrder1::parameter_update_len_call_prob(const vector<int> &lc_count, int c, double gb[2]){
    int idx,idx0,idx1;
    int ell;
    int n=0,nn=0;
    double b0,b1;
    double gw,w;
    vector<double> X, XX;
    vector<double> y, yy;

    // computing the parameters of generalized linear model
    for (int strand=1; strand>=0; strand--){
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            // compute ell
            ell=SemiHomopolymerAlignmentSpace::state_length[pi];

            // pack up X and y
            idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                pi*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                ell;
            idx1=(ell-1)*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                 ell;
            if (lc_count[idx0]>100){
                int dk;
                double tot=0,Z=0;
                double tau=100;
                for (int k=0;k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                        pi*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                        k;
                    dk=abs(k-ell);
                    if (dk<ell){
                        tot+=dk*(lc_count[idx0]+exp(log(tau)-dk*log(tau)));
                        Z+=lc_count[idx0]+exp(log(tau)-dk*log(tau));
                    }else{
                        tot+=dk*lc_count[idx0];
                        Z+=lc_count[idx0];
                    }
                }
                X.push_back(ell);
                y.push_back(tot/Z);
                n++;
            }

            // pack up XX and yy
            int inst=0;
            for (int k=0; k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                     pi*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                     k;
                inst+=lc_count[idx0];
            }
            if (inst>0){
                int dk;
                double tot=0,Z=0;
                for (int k=0;k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                    idx0=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                        pi*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX+
                        k;
                    dk=abs(k-ell);
                    tot+=dk*lc_count[idx0];
                    Z+=lc_count[idx0];
                }
                XX.push_back(ell);
                yy.push_back(tot/Z);
                nn++;
            }

            // compute the parameters of generalized linear model
            if (pi==(int)SemiHomopolymerAlignmentSpace::A20 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::C20 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::G20 ||
                    pi==(int)SemiHomopolymerAlignmentSpace::T20){

                double b[2];
                GLM::irls(X.data(), y.data(), n, b, 100, 1e-7);

                if (pi==(int)SemiHomopolymerAlignmentSpace::A20){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        (int)NucleotideSpace::A;
                }
                if (pi==(int)SemiHomopolymerAlignmentSpace::C20){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        (int)NucleotideSpace::C;
                }
                if (pi==(int)SemiHomopolymerAlignmentSpace::G20){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        (int)NucleotideSpace::G;
                }
                if (pi==(int)SemiHomopolymerAlignmentSpace::T20){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        (int)NucleotideSpace::T;
                }

                w=GLM::sse(XX.data(),yy.data(),nn,b);
                gw=GLM::sse(XX.data(),yy.data(),nn,gb);
                double rho=gw/(w+gw);
                glm_poisson_b0[idx]=rho*b[0]+(1-rho)*gb[0];
                glm_poisson_b1[idx]=rho*b[1]+(1-rho)*gb[1];

                if (glm_poisson_b0[idx]!=glm_poisson_b0[idx]){
                    glm_poisson_b0[idx]=gb[0];
                }
                if (glm_poisson_b1[idx]!=glm_poisson_b1[idx]){
                    glm_poisson_b1[idx]=gb[1];
                }

                X.clear();
                y.clear();
                n=0;
                XX.clear();
                yy.clear();
                nn=0;
            }
        }

        // compute the probability
        double lambda;
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            ell=SemiHomopolymerAlignmentSpace::state_length[pi];

            if (pi>=(int)SemiHomopolymerAlignmentSpace::A1 && pi<=(int)SemiHomopolymerAlignmentSpace::A20){
                idx=strand*cycles*NucleotideSpace::NucleotideSize+
                    c*NucleotideSpace::NucleotideSize+
                    (int)NucleotideSpace::A;
            }
            if (pi>=(int)SemiHomopolymerAlignmentSpace::C1 && pi<=(int)SemiHomopolymerAlignmentSpace::C20){
                idx=strand*cycles*NucleotideSpace::NucleotideSize+
                    c*NucleotideSpace::NucleotideSize+
                    (int)NucleotideSpace::C;
            }
            if (pi>=(int)SemiHomopolymerAlignmentSpace::G1 && pi<=(int)SemiHomopolymerAlignmentSpace::G20){
                idx=strand*cycles*NucleotideSpace::NucleotideSize+
                    c*NucleotideSpace::NucleotideSize+
                    (int)NucleotideSpace::G;
            }
            if (pi>=(int)SemiHomopolymerAlignmentSpace::T1 && pi<=(int)SemiHomopolymerAlignmentSpace::T20){
                idx=strand*cycles*NucleotideSpace::NucleotideSize+
                    c*NucleotideSpace::NucleotideSize+
                    (int)NucleotideSpace::T;
            }
            b0=glm_poisson_b0[idx];
            b1=glm_poisson_b1[idx];
            // poisson mean
            lambda=GLM::invlinkfunc(b0+b1*ell);

            // compute the probability
            for (int k=0; k<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; k++){
                int dk=abs(ell-k);
                double lp=dk*log(lambda)-lambda-lgamma(dk+1);
                idx=length_call_table_index(strand, c, pi, k);
                if (dk==0){
                    log_prob_length_call[idx]=lp;
                }else{
                    log_prob_length_call[idx]=log(0.5)+lp;
                }
            }
        }
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_print
 */
void SemiHomopolymerGHMMOrder1::parameter_print(){
    int idx;
    char buffer[2048];

    // first positive strand, next negative strand
    for (int c=0; c<cycles; c++){
        cout<<"<<<Cycle #"<<c<<">>>"<<endl;
        for (int strand=1; strand>=0; strand--){
            string direction=(strand==1)?"Positive":"Negative";
            cout<<"{{{Strand: "<<direction<<"}}}"<<endl;
            // title of state transition information
            cout<<"[#1: Probability of Hidden State Transition]"<<endl;
            // state transition probability
            sprintf(buffer, "      %10s","SH");
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::I; pi++){
                if (pi<20){
                    string item="A"+std::to_string(pi+1);
                    sprintf(buffer,"%s %10s",buffer,item.c_str());
                }else if (pi<40){
                    string item="C"+std::to_string(pi-20+1);
                    sprintf(buffer,"%s %10s",buffer,item.c_str());
                }else if (pi<60){
                    string item="G"+std::to_string(pi-40+1);
                    sprintf(buffer,"%s %10s",buffer,item.c_str());
                }else if (pi<80){
                    string item="T"+std::to_string(pi-60+1);
                    sprintf(buffer,"%s %10s",buffer,item.c_str());
                }else{
                    sprintf(buffer,"%s %10s",buffer,"I");
                }
            }
            sprintf(buffer,"%s %10s %10s",buffer,"ST","E");
            cout<<buffer<<endl;
            // B
            idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::B);
            sprintf(buffer,"%3s   %10g","B",exp(log_prob_state_trans[idx]));
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::B);
                sprintf(buffer,"%s %10g",buffer,exp(log_prob_state_trans[idx]));
            }
            cout<<buffer<<endl;
            // SH
            idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::SH);
            sprintf(buffer,"%3s   %10g","SH",exp(log_prob_state_trans[idx]));
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::SH);
                sprintf(buffer,"%s %10g",buffer,exp(log_prob_state_trans[idx]));
            }
            cout<<buffer<<endl;
            // A1 to ST
            for (int ppi=(int)SemiHomopolymerAlignmentSpace::A1; ppi<=(int)SemiHomopolymerAlignmentSpace::ST; ppi++){
                if (ppi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                if (ppi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                sprintf(buffer,"%3s   %10g",SemiHomopolymerAlignmentSpace::idx2state[ppi].c_str(),exp(LOGZERO));
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                    if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                    if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                    idx=state_trans_table_index(strand, c, pi, ppi);
                    sprintf(buffer, "%s %10g", buffer, exp(log_prob_state_trans[idx]));
                }
                cout<<buffer<<endl;
            }

            // title of base-calling table
            cout<<"[#2: Base-Call Probability]"<<endl;
            sprintf(buffer,"state pos qual %10s %10s %10s %10s","A","C","G","T");
            cout<<buffer<<endl;
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
                for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                    for (int qual=0; qual<NucleotideSpace::QualityScoreSlot; qual++){
                        sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                        sprintf(buffer,"%s %d %d", buffer, pos, qual);
                        for (int beta=HomopolymerSpace::A; beta<=HomopolymerSpace::T; beta++){
                            idx=base_call_table_index(strand, c, pi, pos, qual, beta);
                            sprintf(buffer,"%s %10g",buffer,exp(log_prob_base_call[idx]));
                        }
                        cout<<buffer<<endl;
                    }
                }
            }

            // title of quality-calling table
            cout<<"[#3: Quality-Score-Call Probability]"<<endl;
            sprintf(buffer,"state pos");
            for (int qual=0; qual<HomopolymerSpace::QualityScoreSlot; qual++){
                sprintf(buffer, "%s %10d", buffer, qual);
            }
            cout<<buffer<<endl;
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
                for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                    sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                    sprintf(buffer, "%s %d", buffer, pos);
                    for (int qual=0; qual<HomopolymerSpace::QualityScoreSlot; qual++){
                        idx=quality_call_table_index(strand, c, pi, pos, qual);
                        sprintf(buffer, "%s %10g",buffer,exp(log_prob_qual_call[idx]));
                    }
                    cout<<buffer<<endl;
                }
            }

            // title of length-calling table
            cout<<"[#4: Length-Call Probability]"<<endl;
            idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::A;
            cout<<"A: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
            idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::C;
            cout<<"C: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
            idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::G;
            cout<<"G: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
            idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::T;
            cout<<"T: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;

            sprintf(buffer,"state");
            for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                sprintf(buffer,"%s %10d",buffer, kappa);
            }
            cout<<buffer<<endl;
            for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
                sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                    idx=length_call_table_index(strand, c , pi, kappa);
                    sprintf(buffer, "%s %10g",buffer,exp(log_prob_length_call[idx]));
                }
                cout<<buffer<<endl;
            }
            cout<<endl<<endl;
        }
    }
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_print
 * @param filename
 */
void SemiHomopolymerGHMMOrder1::parameter_print(string filename){
    ofstream output;
    output.open(filename);

    {
        int idx;
        char buffer[2048];

        // first positive strand, next negative strand
        for (int c=0; c<cycles; c++){
            output<<"<<<Cycle #"<<c<<">>>"<<endl;
            for (int strand=1; strand>=0; strand--){
                string direction=(strand==1)?"Positive":"Negative";
                output<<"{{{Strand: "<<direction<<"}}}"<<endl;
                // title of state transition information
                output<<"[#1: Probability of Hidden State Transition]"<<endl;
                // state transition probability
                sprintf(buffer, "      %10s","SH");
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::I; pi++){
                    if (pi<20){
                        string item="A"+std::to_string(pi+1);
                        sprintf(buffer,"%s %10s",buffer,item.c_str());
                    }else if (pi<40){
                        string item="C"+std::to_string(pi-20+1);
                        sprintf(buffer,"%s %10s",buffer,item.c_str());
                    }else if (pi<60){
                        string item="G"+std::to_string(pi-40+1);
                        sprintf(buffer,"%s %10s",buffer,item.c_str());
                    }else if (pi<80){
                        string item="T"+std::to_string(pi-60+1);
                        sprintf(buffer,"%s %10s",buffer,item.c_str());
                    }else{
                        sprintf(buffer,"%s %10s",buffer,"I");
                    }
                }
                sprintf(buffer,"%s %10s %10s",buffer,"ST","E");
                output<<buffer<<endl;
                // B
                idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::B);
                sprintf(buffer,"%3s   %10g","B",exp(log_prob_state_trans[idx]));
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                    if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                    if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                    idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::B);
                    sprintf(buffer,"%s %10g",buffer,exp(log_prob_state_trans[idx]));
                }
                output<<buffer<<endl;
                // SH
                idx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::SH);
                sprintf(buffer,"%3s   %10g","SH",exp(log_prob_state_trans[idx]));
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                    if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                    if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                    idx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::SH);
                    sprintf(buffer,"%s %10g",buffer,exp(log_prob_state_trans[idx]));
                }
                output<<buffer<<endl;
                // A1 to ST
                for (int ppi=(int)SemiHomopolymerAlignmentSpace::A1; ppi<=(int)SemiHomopolymerAlignmentSpace::ST; ppi++){
                    if (ppi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                    if (ppi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                    sprintf(buffer,"%3s   %10g",SemiHomopolymerAlignmentSpace::idx2state[ppi].c_str(),exp(LOGZERO));
                    for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::E; pi++){
                        if (pi==(int)SemiHomopolymerAlignmentSpace::B) continue;
                        if (pi==(int)SemiHomopolymerAlignmentSpace::SH) continue;
                        idx=state_trans_table_index(strand, c, pi, ppi);
                        sprintf(buffer, "%s %10g", buffer, exp(log_prob_state_trans[idx]));
                    }
                    output<<buffer<<endl;
                }

                // title of base-calling table
                output<<"[#2: Base-Call Probability]"<<endl;
                sprintf(buffer,"state pos qual %10s %10s %10s %10s","A","C","G","T");
                output<<buffer<<endl;
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
                    for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                        for (int qual=0; qual<NucleotideSpace::QualityScoreSlot; qual++){
                            sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                            sprintf(buffer,"%s %d %d", buffer, pos, qual);
                            for (int beta=HomopolymerSpace::A; beta<=HomopolymerSpace::T; beta++){
                                idx=base_call_table_index(strand, c, pi, pos, qual, beta);
                                sprintf(buffer,"%s %10g",buffer,exp(log_prob_base_call[idx]));
                            }
                            output<<buffer<<endl;
                        }
                    }
                }

                // title of quality-calling table
                output<<"[#3: Quality-Score-Call Probability]"<<endl;
                sprintf(buffer,"state pos");
                for (int qual=0; qual<HomopolymerSpace::QualityScoreSlot; qual++){
                    sprintf(buffer, "%s %10d", buffer, qual);
                }
                output<<buffer<<endl;
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::ST; pi++){
                    for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                        sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                        sprintf(buffer, "%s %d", buffer, pos);
                        for (int qual=0; qual<HomopolymerSpace::QualityScoreSlot; qual++){
                            idx=quality_call_table_index(strand, c, pi, pos, qual);
                            sprintf(buffer, "%s %10g",buffer,exp(log_prob_qual_call[idx]));
                        }
                        output<<buffer<<endl;
                    }
                }

                // title of length-calling table
                output<<"[#4: Length-Call Probability]"<<endl;
                idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::A;
                output<<"A: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
                idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::C;
                output<<"C: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
                idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::G;
                output<<"G: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;
                idx=strand*cycles*NucleotideSpace::NucleotideSize+c*NucleotideSpace::NucleotideSize+(int)NucleotideSpace::T;
                output<<"T: "<<glm_poisson_b0[idx]<<" "<<glm_poisson_b1[idx]<<endl;

                sprintf(buffer,"state");
                for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                    sprintf(buffer,"%s %10d",buffer, kappa);
                }
                output<<buffer<<endl;
                for (int pi=0; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
                    sprintf(buffer,"%3s  ",SemiHomopolymerAlignmentSpace::idx2state[pi].c_str());
                    for (int kappa=0; kappa<HomopolymerSpace::HomopolymerSizeMax; kappa++){
                        idx=length_call_table_index(strand, c , pi, kappa);
                        sprintf(buffer, "%s %10g",buffer,exp(log_prob_length_call[idx]));
                    }
                    output<<buffer<<endl;
                }
                output<<endl<<endl;
            }
        }
    }

    output.close();
}

/**
 * @brief SemiHomopolymerGHMMOrder1::parameter_print_binary
 * @param filename
 */
void SemiHomopolymerGHMMOrder1::parameter_print_binary(string filename){
    ofstream output;
    output.open(filename, ios::out | ios::binary);

    {
        int idx;
        double temp;

        output.write(reinterpret_cast<char*>(&cycles), sizeof(cycles));
        // first positive strand, next negative strand
        for (int c=0; c<cycles; c++){
            for (int strand=1; strand>=0; strand--){
                // state transition probability
                int size_log_prob_state_trans=size_prob_state_trans_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                for (int i=0; i<size_log_prob_state_trans; i++){
                    idx=strand*cycles*size_log_prob_state_trans+
                        c*size_log_prob_state_trans+
                        i;
                    temp=exp(log_prob_state_trans[idx]);
                    output.write(reinterpret_cast<char*>(&temp), sizeof(double));
                }

                // title of base-calling table
                int size_log_prob_base_call=size_prob_base_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                for (int i=0; i<size_log_prob_base_call; i++){
                    idx=strand*cycles*size_log_prob_base_call+
                        c*size_log_prob_base_call+
                        i;
                    temp=exp(log_prob_base_call[idx]);
                    output.write(reinterpret_cast<char*>(&temp), sizeof(double));
                }

                // title of quality-calling table
                int size_log_prob_qual_call=size_prob_qual_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                for (int i=0; i<size_log_prob_qual_call; i++){
                    idx=strand*cycles*size_log_prob_qual_call+
                        c*size_log_prob_qual_call+
                        i;
                    temp=exp(log_prob_qual_call[idx]);
                    output.write(reinterpret_cast<char*>(&temp), sizeof(double));
                }

                // title of length-calling table
                for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        i;
                    output.write(reinterpret_cast<char*>(&glm_poisson_b0[idx]), sizeof(double));
                    output.write(reinterpret_cast<char*>(&glm_poisson_b1[idx]), sizeof(double));
                }
                int size_log_prob_len_call=size_prob_length_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                for (int i=0; i<size_log_prob_len_call; i++){
                    idx=strand*cycles*size_log_prob_len_call+
                        c*size_log_prob_len_call+
                        i;
                    temp=exp(log_prob_length_call[idx]);
                    output.write(reinterpret_cast<char*>(&temp), sizeof(double));
                }
            }
        }
    }

    output.close();
}

/**
 * @brief SemiHomopolymerGHMMOrder1::load_model_configure
 * @param filename
 */
void SemiHomopolymerGHMMOrder1::load_model_configure(string filename){
    ifstream input;
    input.open(filename, ios::in | ios::binary);

    {
        input.read(reinterpret_cast<char*>(&cycles), sizeof(cycles));
        // initialize model
        initialization();

        double temp;
        for (int c=0; c<cycles; c++){
            for (int strand=1; strand>=0; strand--){
                int idx;
                int size_log_prob_state_trans=size_prob_state_trans_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                int size_log_prob_base_call=size_prob_base_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                int size_log_prob_qual_call=size_prob_qual_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                int size_log_prob_len_call=size_prob_length_call_table/(1<<HomopolymerSpace::STRANDBITSIZE)/cycles;
                // state transition probability
                for (int i=0; i<size_log_prob_state_trans; i++){
                    idx=strand*cycles*size_log_prob_state_trans+
                        c*size_log_prob_state_trans+
                        i;

                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    if (temp==0){
                        log_prob_state_trans[idx]=LOGZERO;
                    }else{
                        log_prob_state_trans[idx]=log(temp);
                    }
                }
                // base calling probability
                for (int i=0; i<size_log_prob_base_call; i++){
                    idx=strand*cycles*size_log_prob_base_call+
                        c*size_log_prob_base_call+
                        i;

                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    if (temp==0){
                        log_prob_base_call[idx]=LOGZERO;
                    }else{
                        log_prob_base_call[idx]=log(temp);
                    }
                }
                // quality calling probability
                for (int i=0; i<size_log_prob_qual_call; i++){
                    idx=strand*cycles*size_log_prob_qual_call+
                        c*size_log_prob_qual_call+
                        i;

                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    if (temp==0){
                        log_prob_qual_call[idx]=LOGZERO;
                    }else{
                        log_prob_qual_call[idx]=log(temp);
                    }
                }
                // length calling probability
                for (int i=0; i<NucleotideSpace::NucleotideSize; i++){
                    idx=strand*cycles*NucleotideSpace::NucleotideSize+
                        c*NucleotideSpace::NucleotideSize+
                        i;
                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    glm_poisson_b0[idx]=temp;
                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    glm_poisson_b1[idx]=temp;
                }
                for (int i=0; i<size_log_prob_len_call; i++){
                    idx=strand*cycles*size_log_prob_len_call+
                        c*size_log_prob_len_call+
                        i;

                    input.read(reinterpret_cast<char*>(&temp), sizeof(temp));
                    if (temp==0){
                        log_prob_length_call[idx]=LOGZERO;
                    }else{
                        log_prob_length_call[idx]=log(temp);
                    }
                }
            }
        }
    }

    input.close();
}

/**
 * @brief SemiHomopolymerGHMMOrder1::compute_banded_alignment
 * @param target
 * @param query
 * @param quality
 * @param strand
 * @param alignment
 * @param band
 * @return
 */
double SemiHomopolymerGHMMOrder1::compute_banded_alignment(const string &id, const string &target, const string &query,
                                                         const string &quality, int strand,
                                                         SemiHomopolymerAlignment &alignment,
                                                         int band){
    // 1...
    // banding location
    int t0, t1, h0, h1, q0, q1;
    NucleotideAlignmentMethod nucl_util;
    nucl_util.find_exact_match_segment_chain(target, query, t0, t1, q0, q1);
    string clipped_target=target.substr(t0, t1-t0+1);
    // mapping (t0,t1) to (h0,h1)
    HomopolymerSpace::HomopolymerSequence T(clipped_target);
    h0=0;
    h1=T.len-1;

    // length of target homopolymer sequence
    int m=T.len;
    int mm=target.length();
    int dm=t1-t0+1;
    // length of query nucleotide sequence
    int n=query.length();
    int dn=q1-q0+1;
    // query homopolymer sequence
    HomopolymerSpace::HomopolymerSequence Q(query);
    int nh=0;   // maximal homopolymer run
    for (int j=0; j<Q.len; j++){
        if (nh<Q.ell[j]){
            nh=Q.ell[j];
        }
    }


    // 2...
    // alignment band
    int w=(int)(band*n/100.0);
    //int w_half=w;
    int w_half=(int)(w/2.);
    int w_more=(abs(dm-dn)>nh)?abs(dm-dn):nh;
    if (w_half<w_more){
        w_half+=w_more;
    }

    // NOTE:
    // [1] V1 is to restore the Viterbi score V(pi,i,j), the score of the optimal alignment up to
    //     subsequences Target[0..i] and Query[0..j] with the hidden state pi.
    // [2] PPI is to restore the previous hidden state ppi which transites next to the optimal alignment
    //     hidden state pi.
    // [3] EL is to restore the length of emitted sequence at the hidden state pi.

    // Viterbi matrix
//    int V0_SIZE=SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
//                m*n*SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX;
//    double V0[V0_SIZE];
//    memset(V0, LOGZERO, V0_SIZE*sizeof(double));

    int V1_SIZE=SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*m*n;
    double* V1=new double[V1_SIZE];
    for (int i=0; i<V1_SIZE; i++){
        V1[i]=LOGZERO;
    }

    // matrix to restore the information
    int PPI_SIZE=SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*m*n;
    int *PPI=new int[PPI_SIZE];
    for (int i=0; i<PPI_SIZE; i++){
        PPI[i]=-1;
    }

    int EL_SIZE=SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*m*n;
    int *EL=new int[EL_SIZE];
    for (int i=0; i<EL_SIZE; i++){
        EL[i]=-1;
    }

    // NOTE:
    // [1] Q0 and Q1 are to restore the lower and upper positions of the alignment band at each
    //     hidden state.

    // band specification
    vector<int> Q0(m,0);
    vector<int> Q1(m,0);
    // compute the band matrix
    for (int i=h0; i<=h1; i++){
        int qc=int(q0+(T.t1[i]+1)*dn/(dm+0.));
        Q0[i]=qc-w_half;
        Q1[i]=qc+w_half;
        if (Q0[i]<q0) Q0[i]=q0; // origin point is (h0,q0)
        if (Q1[i]>q1) Q1[i]=q1; // ending point is (h1,q1)
    }

    // NOTE:
    // [1] K0 and K1 are to restore the possible length range of emitted sequence at each
    //     hidden state.
    // [2] KK is to restore the possible length deviation from the hidden state.

    // run-length range matrix
    // K0 is lower bound, K1 is upper bound
    vector<int> K0(m*n,0);
    vector<int> K1(m*n,0);
    vector<int> KK(m,0);
    // compute KK
    for (int hi=h0; hi<=h1; hi++){
        string alpha=T.alpha[hi];
        int ell=T.ell[hi];
        string state=alpha+std::to_string(ell);
        int pi=SemiHomopolymerAlignmentSpace::state2idx[state];
        int c=(int)((t0+T.t0[hi])*cycles/(mm+0.1));
        int dk=1;
        for (; dk<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX; dk++){
            int idx1=this->length_call_table_index(strand, c, pi, ell-dk);
            int idx2=this->length_call_table_index(strand, c, pi, ell+dk);
            double lpp=0;
            if (ell-dk>=0)
                lpp+=log_prob_length_call[idx1];
            if (ell+dk<SemiHomopolymerAlignmentSpace::HOMOPOLYMERSIZEMAX)
                lpp+=log_prob_length_call[idx2];
            if (lpp<log(1e-7)) {
                break;
            }
        }
        KK[hi]=dk;
    }
    // compute the range matrix
    for (int hi=h0; hi<=h1; hi++){
        int qc=int(q0+(T.t1[hi]+1)*dn/(dm+0.));
        for (int qj=Q0[hi]; qj<=Q1[hi]; qj++){
            K0[hi*n+qj]=qj-(qc-1+w_half);
            K1[hi*n+qj]=qj-(qc-1-w_half);
            if (K0[hi*n+qj]<T.ell[hi]-KK[hi]) K0[hi*n+qj]=T.ell[hi]-KK[hi];
            if (K0[hi*n+qj]<0) K0[hi*n+qj]=0;
            if (K0[hi*n+qj]>K1[hi*n+qj]) K1[hi*n+qj]=K0[hi*n+qj]+1;
            if (K1[hi*n+qj]>T.ell[hi]+KK[hi]) K1[hi*n+qj]=T.ell[hi]+KK[hi];
            if (K1[hi*n+qj]>qj-q0+1) K1[hi*n+qj]=qj-q0+1;
        }
    }

    // NOTE:
    // [1] To compute the emission probability starting from the Beginning state to
    //     (q0-1)th query position.  The subsequence Query[0..q0-1] is a soft clipped sequence
    //     which thought to be random emission.

    // TODO: head random emission if there is
    for (int j=0; j<q0; j++){
        double viterbi_score=0;
        int ppi=-1;

        // sequence emission
        int bidx=base_call_table_index(strand, 0, (int)SemiHomopolymerAlignmentSpace::SH,
                                       0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize ,NucleotideSpace::nucl2idx[query[j]]);
        viterbi_score+=log_prob_base_call[bidx];

        // quality emission
        int qidx=quality_call_table_index(strand, 0, (int)SemiHomopolymerAlignmentSpace::SH,
                                          0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize);
        viterbi_score+=log_prob_qual_call[qidx];

        // state probability
        if (j==0){
            // hidden state transition
            int pidx=state_trans_table_index(strand, 0, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::B);
            viterbi_score+=log_prob_state_trans[pidx];
            ppi=(int)SemiHomopolymerAlignmentSpace::B;
        }else{
            // hidden state transition
            int pidx=state_trans_table_index(strand, 0, (int)SemiHomopolymerAlignmentSpace::SH, (int)SemiHomopolymerAlignmentSpace::SH);
            viterbi_score+=log_prob_state_trans[pidx];
            ppi=(int)SemiHomopolymerAlignmentSpace::SH;

            int pv1idx=ppi*m*n+
                       h0*n+
                       j-1;
            viterbi_score+=V1[pv1idx];
        }

        // update V1 matrix
        int v1idx=((int)SemiHomopolymerAlignmentSpace::SH)*m*n+
                  h0*n+
                  j;
        V1[v1idx]=viterbi_score;
        PPI[v1idx]=ppi;
        EL[v1idx]=1;
    }

    // NOTE:
    // [1] The alignment is to start from (h0,q0) to (h1,q1)
    // [2] The searching band on Query is specified by Q0 and Q1
    // [3] At each hidden state, the range of emitted sequence is specified by K0 and K1

    // Middle interval: HMM alignment
    for (int i=h0; i<=h1; i++){
        string state=T.alpha[i]+std::to_string(T.ell[i]);
        int pi=SemiHomopolymerAlignmentSpace::state2idx[state];
        int c=(int)((t0+T.t0[i])*cycles/(mm+0.1));
        for (int j=Q0[i]; j<=Q1[i]; j++){
            double max_viterbi_score=LOGZERO;   // maximal-scored alignment ending at i and j with state pi
            double max_ppi=-1;                  // previous state of maximal-scored alignment
            double max_k=-1;                    // emission length of maximal-scored alignment
            double emission_score;
            double state_score;
            // hidden state <alpha,ell>
            for (int k=K0[i*n+j]; k<=K1[i*n+j]; k++){
                emission_score=0;
                // sequence and quality emission
                for (int jj=j-k+1;jj<=j; jj++){
                    int qual=((int)quality[jj]-33);
                    // sequence emission
                    int bi=NucleotideSpace::nucl2idx[query[jj]];
                    int bidx=base_call_table_index(strand, c, pi, jj-(j-k+1), qual/NucleotideSpace::QualityScoreSlotSize, bi);
                    emission_score+=log_prob_base_call[bidx];
                    // quality emission
                    int qidx=quality_call_table_index(strand, c, pi, jj-(j-k+1), qual/NucleotideSpace::QualityScoreSlotSize);
                    emission_score+=log_prob_qual_call[qidx];
                }

                // length emission
                int lidx=length_call_table_index(strand, c, pi, k);
                emission_score+=log_prob_length_call[lidx];


                // state probability
                if (i==h0){
                    if (j-k+1==q0){
                        state_score=0;
                        // if q0=0, then previous state is B
                        if (q0==0){
                            // hidden state transition
                            int pidx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::B);
                            state_score+=log_prob_state_trans[pidx];
                        }else{  // otherwise, the previous state is S
                            // hidden state transition
                            int pidx=state_trans_table_index(strand, c, pi, (int)SemiHomopolymerAlignmentSpace::SH);
                            state_score+=log_prob_state_trans[pidx];
                            // previous score
                            int pv1idx=((int)SemiHomopolymerAlignmentSpace::SH)*m*n+
                                    h0*n+
                                    j-k;
                            state_score+=V1[pv1idx];
                        }
                        // restore the maximum item
                        if (emission_score+state_score>max_viterbi_score){
                            max_viterbi_score=emission_score+state_score;
                            if (q0==0){
                                max_ppi=SemiHomopolymerAlignmentSpace::B;
                            }else{
                                max_ppi=SemiHomopolymerAlignmentSpace::SH;
                            }
                            max_k=k;
                        }
                    }
                }else{
                    if (j-k>=q0){
                        // hidden state transition and previous score
                        for (int ppi=0; ppi<=(int)SemiHomopolymerAlignmentSpace::I; ppi++){
                            state_score=0;
                            // hidden state transtion
                            int pidx=state_trans_table_index(strand, c, pi, ppi);
                            state_score+=log_prob_state_trans[pidx];
                            // previous score
                            int pv1idx=ppi*m*n+
                                    (i-1)*n+
                                    j-k;
                            state_score+=V1[pv1idx];
                            // restore the maximum item
                            if (emission_score+state_score>max_viterbi_score){
                                max_viterbi_score=emission_score+state_score;
                                max_ppi=ppi;
                                max_k=k;
                            }
                        }
                    }
                }
            }
            // update V1 matrix for hidden state pi
            int v1idx=pi*m*n+i*n+j;
            V1[v1idx]=max_viterbi_score;
            PPI[v1idx]=max_ppi;
            EL[v1idx]=max_k;

            // random insertion state
            if (j-1>=q0){
                max_viterbi_score=LOGZERO;
                max_ppi=-1;
                max_k=-1;
                // sequence emission
                int bidx=base_call_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::I,
                                               0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize, NucleotideSpace::nucl2idx[query[j]]);
                emission_score=log_prob_base_call[bidx];
                // quality emission
                int qidx=quality_call_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::I,
                                                  0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize);
                emission_score+=log_prob_qual_call[qidx];
                // state probability
                for (int ppi=0; ppi<=(int)SemiHomopolymerAlignmentSpace::I; ppi++){
                    state_score=0;
                    // hidden state transition
                    int pidx=state_trans_table_index(strand, c, (int)SemiHomopolymerAlignmentSpace::I, ppi);
                    state_score+=log_prob_state_trans[pidx];
                    // previous score
                    int pv1idx=ppi*m*n+
                            i*n+
                            j-1;
                    state_score+=V1[pv1idx];
                    // restore the maximum item
                    if (emission_score+state_score>max_viterbi_score){
                        max_viterbi_score=emission_score+state_score;
                        max_ppi=ppi;
                        max_k=1;
                    }
                }
                // update V1 matrix for random insertion
                v1idx=(int)(SemiHomopolymerAlignmentSpace::I)*m*n+i*n+j;
                V1[v1idx]=max_viterbi_score;
                PPI[v1idx]=max_ppi;
                EL[v1idx]=max_k;
            }
        }

        // debug
        {
            int ss=0;
            for (int j=Q0[i]; j<=Q1[i]; j++){
                if (V1[pi*m*n+i*n+j]>LOGZERO)
                    ss++;
            }
            if (ss==0){
                cout<<id<<endl;
            }
        }
    }

    // NOTE:
    // [1] The subsequence Query[q1+1..n-1] is random emission

    // TODO: tail random emission if there is
    for (int j=q1+1; j<n; j++){
        double max_viterbi_score=LOGZERO;
        double max_ppi=-1;
        double max_k=-1;
        double emission_score=0;
        double state_score=0;

        // sequence emission
        int bidx=base_call_table_index(strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::ST,
                                       0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize, NucleotideSpace::nucl2idx[query[j]]);
        emission_score+=log_prob_base_call[bidx];
        // quality emission
        int qidx=quality_call_table_index(strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::ST,
                                          0, ((int)quality[j]-33)/NucleotideSpace::QualityScoreSlotSize);
        emission_score+=log_prob_qual_call[qidx];
        // state probability
        if (j==q1+1){
            // hidden state transition and previous score
            for (int ppi=0; ppi<=(int)SemiHomopolymerAlignmentSpace::I; ppi++){
                state_score=0;
                // hidden state transition
                int pidx=state_trans_table_index(strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::ST, ppi);
                state_score+=log_prob_state_trans[pidx];
                // previous score
                int pv1idx=ppi*m*n+
                           h1*n+
                           j-1;
                state_score+=V1[pv1idx];
                // restore the maximum item
                if (emission_score+state_score>max_viterbi_score){
                    max_viterbi_score=emission_score+state_score;
                    max_ppi=ppi;
                    max_k=1;
                }
            }
        }else{
            state_score=0;
            // hidden state transition
            int pidx=state_trans_table_index(strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::ST, (int)SemiHomopolymerAlignmentSpace::ST);
            state_score+=log_prob_state_trans[pidx];
            // previous score
            int pv1idx=((int)SemiHomopolymerAlignmentSpace::ST)*m*n+
                       h1*n+
                       j-1;
            state_score+=V1[pv1idx];
            // restore the maximum item
            if (emission_score+state_score>max_viterbi_score){
                max_viterbi_score=emission_score+state_score;
                max_ppi=(int)SemiHomopolymerAlignmentSpace::ST;
                max_k=1;
            }
        }
        // update V1 matrix
        int v1idx=((int)SemiHomopolymerAlignmentSpace::ST)*m*n+
                  h1*n+
                  j;
        V1[v1idx]=max_viterbi_score;
        PPI[v1idx]=max_ppi;
        EL[v1idx]=max_k;
    }

    // Finally run into the Ending state
    // ending state
    double aln_score=LOGZERO;
    double aln_ending_state=-1;
    for (int ppi=0; ppi<=(int)SemiHomopolymerAlignmentSpace::ST; ppi++){
        double score=0;
        // hidden state transition
        int pidx=state_trans_table_index(strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::E, ppi);
        score+=log_prob_state_trans[pidx];
        // previous score
        int pv1idx=ppi*m*n+
                   h1*n+
                   n-1;
        score+=V1[pv1idx];
        // restore the maximum item
        if (score>aln_score){
            aln_score=score;
            aln_ending_state=ppi;
        }
    }


    // 3...
    // Trace back to pack up the alignment
    alignment.align_name=id;
    alignment.raw_target=target;
    alignment.raw_query=query;
    alignment.raw_quality=quality;
    alignment.query_strand=strand;
    alignment.align_target.reset();
    alignment.align_query.reset();
    alignment.align_status.reset();
    // backtrace to restore the alignment
    // 3.1: Pack up the random tail
    for (int j=n-1; j>=q1+1; j--){
        string subseq=query.substr(j,1);
        string subqual=quality.substr(j,1);
        alignment.align_query.insert_head(subseq,subqual);
        alignment.align_target.insert_head("",0);
        alignment.align_status.insert_head(SemiHomopolymerAlignmentSpace::ST);
    }
    // 3.2: Pack up the alignment part
    int aln_pi;
    int aln_i=h1,aln_j=q1;
    if (q1+1==n){
        aln_pi=aln_ending_state;
    }else{
        int pidx=((int)SemiHomopolymerAlignmentSpace::ST)*m*n+
                 h1*n+
                 q1+1;
        aln_pi=PPI[pidx];
    }
    while(aln_pi!=(int)SemiHomopolymerAlignmentSpace::B &&
          aln_pi!=(int)SemiHomopolymerAlignmentSpace::SH){
        // current hidden state is not random insertion
        if (aln_pi<(int)SemiHomopolymerAlignmentSpace::I){
            // target
            alignment.align_target.insert_head(T.alpha[aln_i], T.ell[aln_i]);
            // query
            int elidx=aln_pi*m*n+
                      aln_i*n+
                      aln_j;
            int k=EL[elidx];
            alignment.align_query.insert_head(query.substr(aln_j-k+1,k), quality.substr(aln_j-k+1,k));
            // alignment state
            alignment.align_status.insert_head((SemiHomopolymerAlignmentSpace::AlignmentState)aln_pi);
            // last state and position
            int v1idx=aln_pi*m*n+
                      aln_i*n+
                      aln_j;
            aln_pi=PPI[v1idx];
            aln_i-=1;
            aln_j-=k;
        }
        // current hidden state is random insertion
        if (aln_pi==(int)SemiHomopolymerAlignmentSpace::I){
            // target
            alignment.align_target.insert_head("",0);
            // query
            int elidx=aln_pi*m*n+
                      aln_i*n+
                      aln_j;
            int k=EL[elidx];
            alignment.align_query.insert_head(query.substr(aln_j-k+1,k), quality.substr(aln_j-k+1,k));
            // alignment state
            alignment.align_status.insert_head((SemiHomopolymerAlignmentSpace::AlignmentState)aln_pi);
            // last state and position
            int v1idx=aln_pi*m*n+
                      aln_i*n+
                      aln_j;
            aln_pi=PPI[v1idx];
            aln_j-=k;
        }
    }
    // 3.3: Pack up the random head
    for (int j=q0-1; j>=0; j--){
        string subseq=query.substr(j,1);
        string subqual=quality.substr(j,1);
        alignment.align_query.insert_head(subseq,subqual);
        alignment.align_target.insert_head("",0);
        alignment.align_status.insert_head(SemiHomopolymerAlignmentSpace::SH);
    }

    // clear out memory allocation
    delete [] V1;
    delete [] PPI;
    delete [] EL;

    return aln_score;
}


double SemiHomopolymerGHMMOrder1::compute_banded_realignment(const SemiHomopolymerAlignment &old_align, SemiHomopolymerAlignment &new_align, int band){
    return compute_banded_alignment(old_align.align_name, old_align.raw_target, old_align.raw_query,
                                    old_align.raw_quality, old_align.query_strand,
                                    new_align, band);
}

/**
 * @brief SemiHomopolymerGHMMOrder1::compute_log_likelihood_score
 * @param align
 * @return
 */
double SemiHomopolymerGHMMOrder1::compute_log_likelihood_score(const SemiHomopolymerAlignment &align){
    double lls=0;
    double es=0;    // score of emission part
    double ss=0;    // score of state part
    int pi,ppi;     // pi is current state, ppi is previous state
    string qs,qq;
    int idx;
    int c,mm;
    int t0=0,t1=0;
    int k;

    // length of target sequence
    mm=align.raw_target.size();

    // scan from first position to last position
    for (int i=0; i<align.align_status.len; i++){
        ss=0;
        es=0;
        // update position information
        t0=t1;
        t1+=align.align_target.ell[i];
        c=(int)(t0*cycles/(mm+0.1));

        // query sequence and quality information
        qs=align.align_query.seq[i];
        qq=align.align_query.qual[i];

        // first position is linked from B
        pi=(int)align.align_status.status[i];
        if (i==0){
            // state transition
            ppi=(int)SemiHomopolymerAlignmentSpace::B;
        }else{
            // state transition
            ppi=(int)align.align_status.status[i-1];
        }
        idx=state_trans_table_index(align.query_strand, c, pi, ppi);
        ss+=log_prob_state_trans[idx];

        // base emission
        for (int j=0; j<(int)qs.length(); j++){
            int beta=NucleotideSpace::nucl2idx[qs[j]];
            idx=base_call_table_index(align.query_strand, c, pi, j, ((int)qq[j]-33)/NucleotideSpace::QualityScoreSlotSize, beta);
            es+=log_prob_base_call[idx];
        }
        // quality emission
        for (int j=0; j<(int)qq.length(); j++){
            idx=quality_call_table_index(align.query_strand, c, pi, j, ((int)qq[j]-33)/NucleotideSpace::QualityScoreSlotSize);
            es+=log_prob_qual_call[idx];
        }
        // length emission
        k=qs.length();
        idx=length_call_table_index(align.query_strand, c, pi, k);
        es+=log_prob_length_call[idx];
        lls+=ss+es;
    }

    // run into ending status
    ppi=(int)align.align_status.status[align.align_status.len-1];
    idx=state_trans_table_index(align.query_strand, cycles-1, (int)SemiHomopolymerAlignmentSpace::E, ppi);
    lls+=log_prob_state_trans[idx];

    return lls;
}
