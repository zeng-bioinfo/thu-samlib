#ifndef SEMI_HOMO_GHMM_ORDER1_H
#define SEMI_HOMO_GHMM_ORDER1_H

#include "nucl_align.h"
#include "semi_homo_align.h"


const double LOGZERO=-1e+300;


/**
 * @brief The SemiHomopolymerGHMMOrder1 class
 */
class SemiHomopolymerGHMMOrder1:NucleotideAlignmentMethod
{
    public:
        /**
         * @brief cycles
         * Within individual cycle, HMM use individual parameter setting
         */
        int cycles;

    public:
        /**
         * @brief log_prob_state_trans
         * Hidden state transition probability
         * Dimension: [STRAND][CYCLE][STATE][STATE]
         */
        double *log_prob_state_trans;
        int size_prob_state_trans_table;
        /**
         * @brief log_prob_base_call
         * Hidden state's base emission probability
         * Dimension: [STRAND][CYCLE][STATE][POSITION][QUALITY][BASE]
         */
        double *log_prob_base_call;
        int size_prob_base_call_table;
        /**
         * @brief log_prob_length_call
         * Hidden state's length emission probability
         * Dimension: [STRAND][CYCLE][STATE][LENGTH]
         */
        double *log_prob_length_call;
        int size_prob_length_call_table;
        /**
         * @brief log_prob_qual_call
         * Hidden state's quality emission probability
         * Dimension: [STRAND][CYCLES][STATE][POSITION][QUALITY]
         */
        double *log_prob_qual_call;
        int size_prob_qual_call_table;
        /**
         * @brief glm_poisson_b0
         * Generalized linear model's parameter for poisson model
         * \lambda = b0 + b1*length(STATE)
         * Dimension: [STRAND][CYCLE][BASE]
         */
        double *glm_poisson_b0;
        int size_glm_poisson_b0_table;
        /**
         * @brief glm_poisson_b1
         * Generalized linear model's parameter for poisson model
         * \lambda = b0 + b1*length(STATE)
         * Dimension: [STRAND][CYCLE][BASE]
         */
        double *glm_poisson_b1;
        int size_glm_poisson_b1_table;

    public:
        int count_thresh=10;

    public:
        /**
         * @brief parameter_estimate
         * @param align_pool
         */
        void parameter_estimate(SemiHomopolymerAlignmentPool& align_pool, int band, int maxIter, double thres);
        /**
         * @brief parameter_initialize
         * @param align_pool
         * Initially set the parameters of Generalized HMM by using
         * the given alignments
         */
        void parameter_initialize(SemiHomopolymerAlignmentPool& align_pool);
        /**
         * @brief parameter_update
         * @param statistics
         */
        void parameter_update(vector<SemiHomopolymerAlignmentStat>& statistics);
        /**
         * @brief parameter_update_state_trans_prob
         * @param st_count
         * @param c
         */
        void parameter_update_state_trans_prob(const vector<int>& st_count, int c);
        /**
         * @brief parameter_update_base_call_prob
         * @param bc_count
         * @param c
         */
        void parameter_update_base_call_prob(const vector<int>& bc_count, int c);
        /**
         * @brief parameter_update_qual_call_prob
         * @param qc_count
         * @param c
         */
        void parameter_update_qual_call_prob(const vector<int>& qc_count, int c);
        /**
         * @brief parameter_update_len_call_prob
         * @param lc_count
         * @param c
         */
        void parameter_update_len_call_prob(const vector<int>& lc_count, int c, double gb[2]);
        /**
         * @brief parameter_print
         */
        void parameter_print();
        /**
         * @brief parameter_print
         * @param filename
         */
        void parameter_print(string filename);
        /**
         * @brief parameter_print_binary
         * @param filename
         */
        void parameter_print_binary(string filename);

    public:
        /**
         * @brief initialization
         */
        void initialization();
        /**
         * @brief reset
         */
        void reset();

    public:
        /**
         * @brief compute_alignment
         * @param target
         * @param query
         * @param alignment
         */
        void compute_alignment(const char* target, const char* query, SemiHomopolymerAlignment& alignment);
        /**
         * @brief compute_alignment
         * @param target
         * @param query
         * @param alignment
         */
        void compute_alignment(const string& target, const string& query, SemiHomopolymerAlignment& alignment);
        /**
         * @brief compute_realignment
         * @param raw_align
         * @param new_align
         */
        void compute_realignment(const SemiHomopolymerAlignment& old_align, SemiHomopolymerAlignment& new_align);
        /**
         * @brief compute_realignment
         * @param raw_align
         * @param new_align
         */
        void compute_realignment(const NucleotideAlignment& old_align, SemiHomopolymerAlignment& new_align);

    public:
        /**
         * @brief load_model_configure
         * @param filename
         */
        void load_model_configure(string filename);

    public:
        /**
         * @brief compute_banded_alignment
         * @param target
         * @param query
         * @param alignment
         */
        void compute_banded_alignment(const char* id, const char* target, const char* query, SemiHomopolymerAlignment& alignment, int band);
        /**
         * @brief compute_banded_alignment
         * @param target
         * @param query
         * @param alignment
         */
        double compute_banded_alignment(const string& id, const string& target, const string& query, const string& quality, int strand, SemiHomopolymerAlignment& alignment, int band);
        /**
         * @brief compute_banded_realignment
         * @param old_align
         * @param new_align
         */
        double compute_banded_realignment(const SemiHomopolymerAlignment& old_align, SemiHomopolymerAlignment& new_align, int band);
        /**
         * @brief compute_banded_realignment
         * @param old_align
         * @param new_align
         */
        double compute_banded_realignment(const NucleotideAlignment& old_align, SemiHomopolymerAlignment& new_align, int band);

    public:
        /**
         * @brief compute_likelihood
         * @param align
         */
        double compute_log_likelihood_score(const SemiHomopolymerAlignment& align);

    public:
        /**
         * @brief state_trans_table_index
         * @param strand
         * @param cycle
         * @param pi
         * @param ppi
         * @return
         */
        int state_trans_table_index(int strand, int cycle, int pi, int ppi){
            return strand*cycles*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                   cycle*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                   ppi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                   pi;
        }
        /**
         * @brief base_call_table_index
         * @param strand
         * @param cycle
         * @param pi
         * @param qual
         * @param beta
         * @return
         */
        int base_call_table_index(int strand, int cycle, int pi, int pos, int qual, int beta){
            return strand*cycles*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                   cycle*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                   pi*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                   pos*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                   qual*NucleotideSpace::NucleotideSize+
                   beta;
        }
        /**
         * @brief quality_call_table_index
         * @param strand
         * @param cycle
         * @param pi
         * @param qual
         * @return
         */
        int quality_call_table_index(int strand, int cycle, int pi, int pos, int qual){
            return strand*cycles*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                   cycle*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                   pi*HomopolymerSpace::HomopolymerPosMax*NucleotideSpace::QualityScoreSlot+
                   pos*NucleotideSpace::QualityScoreSlot+
                   qual;
        }
        /**
         * @brief length_call_table_index
         * @param strand
         * @param cycle
         * @param pi
         * @param length
         * @return
         */
        int length_call_table_index(int strand, int cycle, int pi, int length){
            return strand*cycles*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerSizeMax+
                   cycle*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerSizeMax+
                   pi*HomopolymerSpace::HomopolymerSizeMax+
                   length;
        }

    public:
        /**
         * @brief SemiHomopolymerGHMMOrder1
         */
        SemiHomopolymerGHMMOrder1();
        /**
         * @brief SemiHomopolymerGHMMOrder1
         * @param _cycles
         */
        SemiHomopolymerGHMMOrder1(int _cycles);
        /**
         * @brief ~SemiHomopolymerGHMMOrder1
         */
        virtual ~SemiHomopolymerGHMMOrder1();
};

#endif // SEMI_HOMO_GHMM_ORDER1_H
