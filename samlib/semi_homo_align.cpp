#include "semi_homo_align.h"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <algorithm>
#include <map>
#include <boost/assign/list_of.hpp>
#include <boost/lexical_cast.hpp>
using namespace std;


namespace SemiHomopolymerAlignmentSpace{
    std::map<string, int> state2idx=boost::assign::map_list_of
            ("A1", int(A1))("A2", int(A2))("A3", int(A3))("A4", int(A4))("A5", int(A5))("A6", int(A6))("A7", int(A7))("A8", int(A8))("A9", int(A9))("A10", int(A10))
            ("A11", int(A11))("A12", int(A12))("A13", int(A13))("A14", int(A14))("A15", int(A15))("A16", int(A16))("A17", int(A17))("A18", int(A18))("A19", int(A19))("A20", int(A20))
            ("C1", int(C1))("C2", int(C2))("C3", int(C3))("C4", int(C4))("C5", int(C5))("C6", int(C6))("C7", int(C7))("C8", int(C8))("C9", int(C9))("C10", int(C10))
            ("C11", int(C11))("C12", int(C12))("C13", int(C13))("C14", int(C14))("C15", int(C15))("C16", int(C16))("C17", int(C17))("C18", int(C18))("C19", int(C19))("C20", int(C20))
            ("G1", int(G1))("G2", int(G2))("G3", int(G3))("G4", int(G4))("G5", int(G5))("G6", int(G6))("G7", int(G7))("G8", int(G8))("G9", int(G9))("G10", int(G10))
            ("G11", int(G11))("G12", int(G12))("G13", int(G13))("G14", int(G14))("G15", int(G15))("G16", int(G16))("G17", int(G17))("G18", int(G18))("G19", int(G19))("G20", int(G20))
            ("T1", int(T1))("T2", int(T2))("T3", int(T3))("T4", int(T4))("T5", int(T5))("T6", int(T6))("T7", int(T7))("T8", int(T8))("T9", int(T9))("T10", int(T10))
            ("T11", int(T11))("T12", int(T12))("T13", int(T13))("T14", int(T14))("T15", int(T15))("T16", int(T16))("T17", int(T17))("T18", int(T18))("T19", int(T19))("T20", int(T20))
            ("I", int(I))("SH", int(SH))("ST", int(ST))("B", int(B))("E", int(E));
    std::map<int, string> idx2state=boost::assign::map_list_of
            (int(A1), "A1")(int(A2), "A2")(int(A3), "A3")(int(A4), "A4")(int(A5), "A5")(int(A6), "A6")(int(A7), "A7")(int(A8), "A8")(int(A9), "A9")(int(A10), "A10")
            (int(A11), "A11")(int(A12), "A12")(int(A13), "A13")(int(A14), "A14")(int(A15), "A15")(int(A16), "A16")(int(A17), "A17")(int(A18), "A18")(int(A19), "A19")(int(A20), "A20")
            (int(C1), "C1")(int(C2), "C2")(int(C3), "C3")(int(C4), "C4")(int(C5), "C5")(int(C6), "C6")(int(C7), "C7")(int(C8), "C8")(int(C9), "C9")(int(C10), "C10")
            (int(C11), "C11")(int(C12), "C12")(int(C13), "C13")(int(C14), "C14")(int(C15), "C15")(int(C16), "C16")(int(C17), "C17")(int(C18), "C18")(int(C19), "C19")(int(C20), "C20")
            (int(G1), "G1")(int(G2), "G2")(int(G3), "G3")(int(G4), "G4")(int(G5), "G5")(int(G6), "G6")(int(G7), "G7")(int(G8), "G8")(int(G9), "G9")(int(G10), "G10")
            (int(G11), "G11")(int(G12), "G12")(int(G13), "G13")(int(G14), "G14")(int(G15), "G15")(int(G16), "G16")(int(G17), "G17")(int(G18), "G18")(int(G19), "G19")(int(G20), "G20")
            (int(T1), "T1")(int(T2), "T2")(int(T3), "T3")(int(T4), "T4")(int(T5), "T5")(int(T6), "T6")(int(T7), "T7")(int(T8), "T8")(int(T9), "T9")(int(T10), "T10")
            (int(T11), "T11")(int(T12), "T12")(int(T13), "T13")(int(T14), "T14")(int(T15), "T15")(int(T16), "T16")(int(T17), "T17")(int(T18), "T18")(int(T19), "T19")(int(T20), "T20")
            (int(I), "I")(int(SH), "SH")(int(ST), "ST")(int(B), "B")(int(E), "E");
    std::map<int, int> state_length=boost::assign::map_list_of
            (int(A1), 1)(int(A2), 2)(int(A3), 3)(int(A4), 4)(int(A5), 5)(int(A6), 6)(int(A7), 7)(int(A8), 8)(int(A9), 9)(int(A10), 10)
            (int(A11), 11)(int(A12), 12)(int(A13), 13)(int(A14), 14)(int(A15), 15)(int(A16), 16)(int(A17), 17)(int(A18), 18)(int(A19), 19)(int(A20), 20)
            (int(C1), 1)(int(C2), 2)(int(C3), 3)(int(C4), 4)(int(C5), 5)(int(C6), 6)(int(C7), 7)(int(C8), 8)(int(C9), 9)(int(C10), 10)
            (int(C11), 11)(int(C12), 12)(int(C13), 13)(int(C14), 14)(int(C15), 15)(int(C16), 16)(int(C17), 17)(int(C18), 18)(int(C19), 19)(int(C20), 20)
            (int(G1), 1)(int(G2), 2)(int(G3), 3)(int(G4), 4)(int(G5), 5)(int(G6), 6)(int(G7), 7)(int(G8), 8)(int(G9), 9)(int(G10), 10)
            (int(G11), 11)(int(G12), 12)(int(G13), 13)(int(G14), 14)(int(G15), 15)(int(G16), 16)(int(G17), 17)(int(G18), 18)(int(G19), 19)(int(G20), 20)
            (int(T1), 1)(int(T2), 2)(int(T3), 3)(int(T4), 4)(int(T5), 5)(int(T6), 6)(int(T7), 7)(int(T8), 8)(int(T9), 9)(int(T10), 10)
            (int(T11), 11)(int(T12), 12)(int(T13), 13)(int(T14), 14)(int(T15), 15)(int(T16), 16)(int(T17), 17)(int(T18), 18)(int(T19), 19)(int(T20), 20)
            (int(I), 0)(int(SH), 0)(int(ST), 0)(int(B), 0)(int(E), 0);
    std::map<int, char> state_alpha=boost::assign::map_list_of
            (int(A1), 'A')(int(A2), 'A')(int(A3), 'A')(int(A4), 'A')(int(A5), 'A')(int(A6), 'A')(int(A7), 'A')(int(A8), 'A')(int(A9), 'A')(int(A10), 'A')
            (int(A11), 'A')(int(A12), 'A')(int(A13), 'A')(int(A14), 'A')(int(A15), 'A')(int(A16), 'A')(int(A17), 'A')(int(A18), 'A')(int(A19), 'A')(int(A20), 'A')
            (int(C1), 'C')(int(C2), 'C')(int(C3), 'C')(int(C4), 'C')(int(C5), 'C')(int(C6), 'C')(int(C7), 'C')(int(C8), 'C')(int(C9), 'C')(int(C10), 'C')
            (int(C11), 'C')(int(C12), 'C')(int(C13), 'C')(int(C14), 'C')(int(C15), 'C')(int(C16), 'C')(int(C17), 'C')(int(C18), 'C')(int(C19), 'C')(int(C20), 'C')
            (int(G1), 'G')(int(G2), 'G')(int(G3), 'G')(int(G4), 'G')(int(G5), 'G')(int(G6), 'G')(int(G7), 'G')(int(G8), 'G')(int(G9), 'G')(int(G10), 'G')
            (int(G11), 'G')(int(G12), 'G')(int(G13), 'G')(int(G14), 'G')(int(G15), 'G')(int(G16), 'G')(int(G17), 'G')(int(G18), 'G')(int(G19), 'G')(int(G20), 'G')
            (int(T1), 'T')(int(T2), 'T')(int(T3), 'T')(int(T4), 'T')(int(T5), 'T')(int(T6), 'T')(int(T7), 'T')(int(T8), 'T')(int(T9), 'T')(int(T10), 'T')
            (int(T11), 'T')(int(T12), 'T')(int(T13), 'T')(int(T14), 'T')(int(T15), 'T')(int(T16), 'T')(int(T17), 'T')(int(T18), 'T')(int(T19), 'T')(int(T20), 'T');
}

/**
 * @brief SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat
 */
SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat(){
    // base call table
    this->base_call_table.resize((1<<HomopolymerSpace::STRANDBITSIZE)*
                                      SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
                                      NucleotideSpace::NucleotideSize, 0);
    // quality score call table
    this->qual_call_table.resize((1<<HomopolymerSpace::STRANDBITSIZE)*
                                      SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
                                      NucleotideSpace::QualityScoreMax, 0);
    // length call table
    this->len_call_table.resize((1<<HomopolymerSpace::STRANDBITSIZE)*
                                     SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
                                     HomopolymerSpace::HomopolymerSizeMax, 0);
    // positional quality-score-related base-calling table
    this->pos_qual_base_call_table.resize((1<<HomopolymerSpace::STRANDBITSIZE)*
                                          SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
                                          HomopolymerSpace::HomopolymerPosMax*
                                          HomopolymerSpace::QualityScoreSlot*
                                          NucleotideSpace::NucleotideSize);
    // state transition table
    this->state_trans_table.resize((1<<HomopolymerSpace::STRANDBITSIZE)*
                                        SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*
                                        SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE, 0);
}

/**
 * @brief SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat
 * @param another
 */
SemiHomopolymerAlignmentStat::SemiHomopolymerAlignmentStat(const SemiHomopolymerAlignmentStat &another){
    this->base_call_table.assign(another.base_call_table.begin(), another.base_call_table.end());
    this->qual_call_table.assign(another.qual_call_table.begin(), another.qual_call_table.end());
    this->len_call_table.assign(another.len_call_table.begin(), another.len_call_table.end());
    this->pos_qual_base_call_table.assign(another.pos_qual_base_call_table.begin(), another.pos_qual_base_call_table.end());
    this->state_trans_table.assign(another.state_trans_table.begin(), another.state_trans_table.end());
}

/**
 * @brief SemiHomopolymerAlignmentStat::~SemiHomopolymerAlignmentStat
 */
SemiHomopolymerAlignmentStat::~SemiHomopolymerAlignmentStat(){
    vector<int>().swap(this->base_call_table);
    vector<int>().swap(this->qual_call_table);
    vector<int>().swap(this->len_call_table);
    vector<int>().swap(this->pos_qual_base_call_table);
    vector<int>().swap(this->state_trans_table);
}

/**
 * @brief SemiHomopolymerAlignmentStat::operator =
 * @param another
 * @return
 */
SemiHomopolymerAlignmentStat& SemiHomopolymerAlignmentStat::operator =(const SemiHomopolymerAlignmentStat &another){
    this->base_call_table.assign(another.base_call_table.begin(), another.base_call_table.end());
    this->qual_call_table.assign(another.qual_call_table.begin(), another.qual_call_table.end());
    this->len_call_table.assign(another.len_call_table.begin(), another.len_call_table.end());
    this->pos_qual_base_call_table.assign(another.pos_qual_base_call_table.begin(), another.pos_qual_base_call_table.end());
    this->state_trans_table.assign(another.state_trans_table.begin(), another.state_trans_table.end());
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
    // update positional quality-score-related base call table
    for (int i=0; i<(int)pos_qual_base_call_table.size(); i++){
        pos_qual_base_call_table[i]+=another.pos_qual_base_call_table[i];
    }

    // update state transition table
    for (int i=0; i<(int)state_trans_table.size(); i++){
        state_trans_table[i]+=another.state_trans_table[i];
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE)) |
                 (ell<<HomopolymerSpace::BETABITSIZE) | (HomopolymerSpace::nucl2idx[beta]);
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE)) |
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE)) |
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::BETABITSIZE)) |
                 (ell<<HomopolymerSpace::BETABITSIZE) | (HomopolymerSpace::nucl2idx[beta]);
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::QUALITYBITSIZE)) |
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
                  HomopolymerSpace::nucl2idx[alpha]<<(HomopolymerSpace::LENGTHBITSIZE+HomopolymerSpace::LENGTHBITSIZE)) |
                 (ell<<HomopolymerSpace::LENGTHBITSIZE) | kappa;
    return len_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::bc_incr1
 * @param strand
 * @param pi
 * @param beta
 */
void SemiHomopolymerAlignmentStat::bc_incr1(int strand, int pi, char beta){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*NucleotideSpace::NucleotideSize+
                 pi*NucleotideSpace::NucleotideSize+
                 NucleotideSpace::nucl2idx[beta];
    base_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::qc_incr1
 * @param strand
 * @param pi
 * @param qual
 */
void SemiHomopolymerAlignmentStat::qc_incr1(int strand, int pi, int qual){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*NucleotideSpace::QualityScoreMax+
                 pi*NucleotideSpace::QualityScoreMax+
                 qual;
    qual_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::pqbc_incr1
 * @param strand
 * @param pos
 * @param qual
 * @param alpha
 * @param beta
 */
void SemiHomopolymerAlignmentStat::pqbc_incr1(int strand, int pi, int pos, int qual, int beta){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 pi*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 pos*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 qual*NucleotideSpace::NucleotideSize+
                 beta;
    pos_qual_base_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::lc_incr1
 * @param strand
 * @param pi
 * @param kappa
 */
void SemiHomopolymerAlignmentStat::lc_incr1(int strand, int pi, int kappa){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerSizeMax+
                 pi*HomopolymerSpace::HomopolymerSizeMax+
                 kappa;
    len_call_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::bc_elem1
 * @param strand
 * @param pi
 * @param beta
 * @return
 */
int SemiHomopolymerAlignmentStat::bc_elem1(int strand, int pi, char beta){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*NucleotideSpace::NucleotideSize+
                 pi*NucleotideSpace::NucleotideSize+
                 NucleotideSpace::nucl2idx[beta];
    return base_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::qc_elem1
 * @param strand
 * @param pi
 * @param qual
 * @return
 */
int SemiHomopolymerAlignmentStat::qc_elem1(int strand, int pi, int qual){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*NucleotideSpace::QualityScoreMax+
                 pi*NucleotideSpace::QualityScoreMax+
                 qual;
    return qual_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::pqbc_elem1
 * @param strand
 * @param pos
 * @param qual
 * @param alpha
 * @param beta
 * @return
 */
int SemiHomopolymerAlignmentStat::pqbc_elem1(int strand, int pi, int pos, int qual, int beta){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 pi*HomopolymerSpace::HomopolymerPosMax*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 pos*HomopolymerSpace::QualityScoreSlot*NucleotideSpace::NucleotideSize+
                 qual*NucleotideSpace::NucleotideSize+
                 beta;
    return pos_qual_base_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::lc_elem1
 * @param strand
 * @param pi
 * @param kappa
 * @return
 */
int SemiHomopolymerAlignmentStat::lc_elem1(int strand, int pi, int kappa){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*HomopolymerSpace::HomopolymerSizeMax+
                 pi*HomopolymerSpace::HomopolymerSizeMax+
                 kappa;
    return len_call_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::st_incr1
 * @param pi
 * @param ppi
 */
void SemiHomopolymerAlignmentStat::st_incr1(int strand, int pi, int ppi){
    uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                 pi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                 ppi;
    state_trans_table[idx]+=1;
}

/**
 * @brief SemiHomopolymerAlignmentStat::st_elem1
 * @param pi
 * @param ppi
 * @return
 */
int SemiHomopolymerAlignmentStat::st_elem1(int strand, int pi, int ppi){
   uint32_t idx=strand*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                pi*SemiHomopolymerAlignmentSpace::ALIGNMENTSTATSIZE+
                ppi;
   return state_trans_table[idx];
}

/**
 * @brief SemiHomopolymerAlignmentStat::print
 * @param c
 */
void SemiHomopolymerAlignmentStat::print(int c){
    // print out the statistics
    char buffer[8192];
    if (c==0){
        // print header
        sprintf(buffer, "cycle strand state pos qual match mismatch");
        for (int alpha=(int)NucleotideSpace::A; alpha<=(int)NucleotideSpace::T; alpha++){
            for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                sprintf(buffer, "%s %c%c", buffer,
                        NucleotideSpace::idx2nucl[alpha],
                        NucleotideSpace::idx2nucl[beta]);
            }
        }
        cout<<buffer<<endl;
    }

    for (int strand=1; strand>=0; strand--){
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual+=HomopolymerSpace::QualityScoreSlotSize){
                    sprintf(buffer, "%d %d %s %d %d", c, strand, idx2state[pi].c_str(), pos+1, qual/HomopolymerSpace::QualityScoreSlotSize);
                    // match and mismatch
                    int match=0;
                    int mismatch=0;
                    for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                        if (state_alpha[pi]!=NucleotideSpace::idx2nucl[beta]){
                            mismatch+=pqbc_elem1(strand, pi, pos,
                                                 qual/HomopolymerSpace::QualityScoreSlotSize,
                                                 beta);
                        }else{
                            match+=pqbc_elem1(strand, pi, pos,
                                              qual/HomopolymerSpace::QualityScoreSlotSize,
                                              beta);
                        }
                    }
                    sprintf(buffer, "%s %d %d", buffer, match, mismatch);
                    // alpha->beta
                    for (int alpha=(int)NucleotideSpace::A; alpha<=(int)NucleotideSpace::T; alpha++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            int count=0;

                            if(NucleotideSpace::idx2nucl[alpha]==state_alpha[pi]){
                                count+=pqbc_elem1(strand, pi, pos,
                                                  qual/HomopolymerSpace::QualityScoreSlotSize,
                                                  beta);
                            }

                            sprintf(buffer, "%s %d", buffer, count);
                        }
                    }

                    cout<<buffer<<endl;
                }
            }
        }
    }
}

/**
 * @brief SemiHomopolymerAlignmentStat::print
 * @param filename
 */
void SemiHomopolymerAlignmentStat::print(ofstream &ofs, int c){
    // print out the statistics
    char buffer[8192];
    if (c==0){
        // print header
        sprintf(buffer, "cycle strand state pos qual match mismatch");
        for (int alpha=(int)NucleotideSpace::A; alpha<=(int)NucleotideSpace::T; alpha++){
            for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                sprintf(buffer, "%s %c%c", buffer,
                        NucleotideSpace::idx2nucl[alpha],
                        NucleotideSpace::idx2nucl[beta]);
            }
        }
        ofs<<buffer<<endl;
    }

    for (int strand=1; strand>=0; strand--){
        for (int pi=(int)SemiHomopolymerAlignmentSpace::A1; pi<=(int)SemiHomopolymerAlignmentSpace::T20; pi++){
            for (int pos=0; pos<HomopolymerSpace::HomopolymerPosMax; pos++){
                for (int qual=0; qual<HomopolymerSpace::QualityScoreMax; qual+=HomopolymerSpace::QualityScoreSlotSize){
                    sprintf(buffer, "%d %d %s %d %d", c, strand, idx2state[pi].c_str(), pos+1, qual/HomopolymerSpace::QualityScoreSlotSize);
                    // match and mismatch
                    int match=0;
                    int mismatch=0;
                    for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                        if (state_alpha[pi]!=NucleotideSpace::idx2nucl[beta]){
                            mismatch+=pqbc_elem1(strand, pi, pos,
                                                 qual/HomopolymerSpace::QualityScoreSlotSize,
                                                 beta);
                        }else{
                            match+=pqbc_elem1(strand, pi, pos,
                                              qual/HomopolymerSpace::QualityScoreSlotSize,
                                              beta);
                        }
                    }
                    sprintf(buffer, "%s %d %d", buffer, match, mismatch);
                    // alpha->beta
                    for (int alpha=(int)NucleotideSpace::A; alpha<=(int)NucleotideSpace::T; alpha++){
                        for (int beta=(int)NucleotideSpace::A; beta<=(int)NucleotideSpace::T; beta++){
                            int count=0;

                            if(NucleotideSpace::idx2nucl[alpha]==state_alpha[pi]){
                                count+=pqbc_elem1(strand, pi, pos,
                                                  qual/HomopolymerSpace::QualityScoreSlotSize,
                                                  beta);
                            }

                            sprintf(buffer, "%s %d", buffer, count);
                        }
                    }

                    ofs<<buffer<<endl;
                }
            }
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

    this->raw_query=homo_align.raw_query;
    this->raw_quality=homo_align.raw_quality;
    this->raw_target=homo_align.raw_target;
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

    this->raw_query=nucl_align.raw_query;
    this->raw_quality=nucl_align.raw_quality;
    this->raw_target=nucl_align.raw_target;

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

    this->raw_query=homo_align.raw_query;
    this->raw_quality=homo_align.raw_quality;
    this->raw_target=homo_align.raw_target;

    return *this;
}

SemiHomopolymerAlignment& SemiHomopolymerAlignment::operator =(const NucleotideAlignment& nucl_align){
    this->align_name=nucl_align.align_name;
    this->query_seq=nucl_align.query_seq;
    this->query_qual=nucl_align.query_qual;
    this->query_strand=nucl_align.query_strand;
    this->target_seq=nucl_align.target_seq;

    this->raw_query=nucl_align.raw_query;
    this->raw_quality=nucl_align.raw_quality;
    this->raw_target=nucl_align.raw_target;

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

    string ss="";          // alignment state in nucleotide space (str)
    SemiHomopolymerAlignmentSpace::AlignmentState as;  // alignment state in nucleotide space

    bool at_head=true;

    // initialization
    qis=""; qiq="";
    tia=""; til=0;
    // scan over the nucleotide alignment
    this->align_name=nucl_align.align_name;
    for (int i=0; i<(int)nucl_align.align_status.length(); i++){
        // update alignment state if it is empty
        if (ss.empty()){
            ss=nucl_align.align_status[i];
        }else if (nucl_align.align_status[i]=='S' ||
                  nucl_align.align_status[i]=='s'){
            ss="S";
        }
        if (nucl_align.align_status[i]!='S' &&
                nucl_align.align_status[i]!='s'){
            at_head=false;
        }
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
            // pack up the alignment state information and complete homopolymer information
            if (ss=="S"){
                // alignment state
                if (at_head){
                    as=SemiHomopolymerAlignmentSpace::SH;
                }else{
                    as=SemiHomopolymerAlignmentSpace::ST;
                }
                // pack up
                for (int j=0; j<(int)qis.size(); j++){
                    this->align_query.push_back(qis.substr(j,1),qiq.substr(j,1));
                    this->align_target.push_back("",0);
                    this->align_status.push_back(as);
                }
            }else{
                string state=tia+std::to_string(til);
                // alignment state
                if (til==0){
                    as=SemiHomopolymerAlignmentSpace::I;
                }else{
                    as=(SemiHomopolymerAlignmentSpace::AlignmentState)SemiHomopolymerAlignmentSpace::state2idx[state];
                }
                // pack up
                this->align_query.push_back(qis,qiq);
                this->align_target.push_back(tia,til);
                this->align_status.push_back(as);
            }

            // reset temporaty buffer
            qis=""; qiq="";
            tia=""; til=0;
            ss="";
        }
    }
}

/**
 * @brief SemiHomopolymerAlignment::print
 */
void SemiHomopolymerAlignment::print(){

    // temporary buffer
    string buffer;

    // print out message
    cout<<this->align_name<<endl;

    // print out query sequence
    for (int i=0; i<align_query.len-1; i++){
        buffer=align_query.seq[i];
        if (align_target.ell[i]>(int)buffer.length()){
            for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                buffer+="-";
            }
        }
        cout<<buffer<<flush;
    }
    buffer=align_query.seq[align_query.len-1];
    if (align_target.ell[align_target.len-1]>(int)buffer.length()){
        for (int j=align_target.ell[align_query.len-1]-buffer.length(); j>0; j--){
            buffer+="-";
        }
    }
    cout<<buffer<<endl;

    // print out query quality score
    for (int i=0; i<align_query.len-1; i++){
        buffer=align_query.qual[i];
        if (align_target.ell[i]>(int)buffer.length()){
            for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                buffer+=char(33);
            }
        }
        cout<<buffer<<flush;
    }
    buffer=align_query.qual[align_query.len-1];
    if (align_target.ell[align_target.len-1]>(int)buffer.length()){
        for (int j=align_target.ell[align_target.len-1]-buffer.length(); j>0; j--){
            buffer+=char(33);
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
                buffer+="-";
            }
        }
        cout<<buffer<<flush;
    }
    buffer="";
    for (int j=0; j<align_target.ell[align_target.len-1]; j++){
        buffer+=align_target.alpha[align_target.len-1];
    }
    if (align_target.ell[align_target.len-1]<(int)align_query.seq[align_query.len-1].length()){
        for (int j=align_query.seq[align_query.len-1].length()-align_target.ell[align_target.len-1]; j>0; j--){
            buffer+="-";
        }
    }
    cout<<buffer<<endl;

    // print out alignment status

    for (int i=0; i<align_status.len; i++){
        buffer="";
        if (align_status.status[i]==SemiHomopolymerAlignmentSpace::SH){
            buffer+="S";
        }else if (align_status.status[i]==SemiHomopolymerAlignmentSpace::ST){
            buffer+="S";
        }else if (align_status.status[i]==SemiHomopolymerAlignmentSpace::I){
            buffer+="I";
        }else{
            string alpha=align_target.alpha[i];
            int ell=align_target.ell[i];
            string query_subseq=align_query.seq[i];
            int k=(ell<query_subseq.length())?ell:query_subseq.length();
            for (int j=0; j<k; j++){
                if (alpha==query_subseq.substr(j,1)){
                    buffer+="M";
                }else{
                    buffer+="X";
                }
            }
            for (int j=k; j<ell; j++){
                buffer+="D";
            }
            for (int j=k; j<(int)query_subseq.length(); j++){
                buffer+="I";
            }
        }
        cout<<buffer<<flush;
    }
    cout<<endl;
}

void SemiHomopolymerAlignment::print(string filename){
    ofstream output;
    output.open(filename, ios::out | ios::app);

    {
        // temporary buffer
        string buffer;

        // print out message
        output<<this->align_name<<endl;

        // print out query sequence
        for (int i=0; i<align_query.len-1; i++){
            buffer=align_query.seq[i];
            if (align_target.ell[i]>(int)buffer.length()){
                for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                    buffer+="-";
                }
            }
            output<<buffer<<flush;
        }
        buffer=align_query.seq[align_query.len-1];
        if (align_target.ell[align_target.len-1]>(int)buffer.length()){
            for (int j=align_target.ell[align_query.len-1]-buffer.length(); j>0; j--){
                buffer+="-";
            }
        }
        output<<buffer<<endl;

        // print out query quality score
        for (int i=0; i<align_query.len-1; i++){
            buffer=align_query.qual[i];
            if (align_target.ell[i]>(int)buffer.length()){
                for (int j=align_target.ell[i]-buffer.length(); j>0; j--){
                    buffer+=char(33);
                }
            }
            output<<buffer<<flush;
        }
        buffer=align_query.qual[align_query.len-1];
        if (align_target.ell[align_target.len-1]>(int)buffer.length()){
            for (int j=align_target.ell[align_target.len-1]-buffer.length(); j>0; j--){
                buffer+=char(33);
            }
        }
        output<<buffer<<endl;

        // print out target sequence
        for (int i=0; i<align_target.len-1; i++){
            buffer="";
            for (int j=0; j<align_target.ell[i]; j++){
                buffer+=align_target.alpha[i];
            }
            if (align_target.ell[i]<(int)align_query.seq[i].length()){
                for (int j=align_query.seq[i].length()-align_target.ell[i]; j>0; j--){
                    buffer+="-";
                }
            }
            output<<buffer<<flush;
        }
        buffer="";
        for (int j=0; j<align_target.ell[align_target.len-1]; j++){
            buffer+=align_target.alpha[align_target.len-1];
        }
        if (align_target.ell[align_target.len-1]<(int)align_query.seq[align_query.len-1].length()){
            for (int j=align_query.seq[align_query.len-1].length()-align_target.ell[align_target.len-1]; j>0; j--){
                buffer+="-";
            }
        }
        output<<buffer<<endl;

        // print out alignment status

        for (int i=0; i<align_status.len; i++){
            buffer="";
            if (align_status.status[i]==SemiHomopolymerAlignmentSpace::SH){
                buffer+="S";
            }else if (align_status.status[i]==SemiHomopolymerAlignmentSpace::ST){
                buffer+="S";
            }else if (align_status.status[i]==SemiHomopolymerAlignmentSpace::I){
                buffer+="I";
            }else{
                string alpha=align_target.alpha[i];
                int ell=align_target.ell[i];
                string query_subseq=align_query.seq[i];
                int k=(ell<query_subseq.length())?ell:query_subseq.length();
                for (int j=0; j<k; j++){
                    if (alpha==query_subseq.substr(j,1)){
                        buffer+="M";
                    }else{
                        buffer+="X";
                    }
                }
                for (int j=k; j<ell; j++){
                    buffer+="D";
                }
                for (int j=k; j<(int)query_subseq.length(); j++){
                    buffer+="I";
                }
            }
            output<<buffer<<flush;
        }
        output<<endl;
    }

    output.close();
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
    vector<SemiHomopolymerAlignmentStat>().swap(stat);
    stat.resize(cycles);

    // target length without spaces
    int len=0;
    for (int i=0; i<align_target.len; i++){
        len+=align_target.ell[i];
    }

    // scan through the alignment
    int t0=0,t1=0;
    for (int i=0; i<align_target.len; i++){
        t0=t1;
        t1+=align_target.ell[i];
        // cycle no.
        int c=(int)(t0*cycles/(len+0.1));

        // alignment state information
        int pi,ppi;
        pi=(int)align_status.status[i];
        if (i==0){
            ppi=(int)SemiHomopolymerAlignmentSpace::AlignmentState::B;
        }else{
            ppi=(int)align_status.status[i-1];
        }
        if (pi<0 || ppi<0) continue;
        stat[c].st_incr1(query_strand, pi, ppi);

        // sequence information
        // symbols and size
        char beta;
        int kappa,qual;

        qual=0;
        kappa=align_query.seq[i].size();
        for (int k=0; k<kappa; k++){
            beta=align_query.seq[i][k];
            qual+=(int)align_query.qual[i][k]-33;
            if (HomopolymerSpace::nucl2idx[beta]<HomopolymerSpace::NucleotideSize){
                // update base-call table
                stat[c].bc_incr1(query_strand, pi, beta);
            }
        }
        // update quality-score-call table
        if (kappa>0)
            stat[c].qc_incr1(query_strand, pi, int(qual/(kappa+0.)));

        // update pos_qual_base_call table
        for (int k=0; k<kappa; k++){
            beta=align_query.seq[i][k];
            qual=(int)align_query.qual[i][k]-33;
            if (NucleotideSpace::nucl2idx[beta]<NucleotideSpace::NucleotideSize){
                stat[c].pqbc_incr1(query_strand, pi, k,
                                   qual/HomopolymerSpace::QualityScoreSlotSize,
                                   NucleotideSpace::nucl2idx[beta]);
            }
        }

        // update length-call table
        stat[c].lc_incr1(query_strand, pi, kappa);
    }

    // run into ending state
    stat[cycles-1].st_incr1(query_strand, (int)SemiHomopolymerAlignmentSpace::AlignmentState::E, (int)align_status.status[align_status.len-1]);
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
    cout<<"Load data from "<<filename<<endl;
    // clear everything
    vector<SemiHomopolymerAlignment>().swap(this->align_pool);

    // first, load in nucleotide space
    NucleotideAlignmentPool nucl_space;
    nucl_space.open(filename);
    // second, covert to homopolymer space
    cout<<"Use "<<nucl_space.align_pool.size()<<" alignments"<<endl;

    for (int i=0; i<(int)nucl_space.align_pool.size(); i++){
        int other=0;
        for (int j=0; j<nucl_space.align_pool[i].raw_target.length(); j++){
            if (nucl_space.align_pool[i].raw_target[j]!='A' &&
                    nucl_space.align_pool[i].raw_target[j]!='C' &&
                    nucl_space.align_pool[i].raw_target[j]!='G' &&
                    nucl_space.align_pool[i].raw_target[j]!='T')
                other++;
        }

        if (other>0) continue;

        SemiHomopolymerAlignment aln1(nucl_space.align_pool[i]);
        this->align_pool.push_back(aln1);
    }
}

/**
 * @brief SemiHomopolymerAlignmentPool::open
 * @param filename
 * @param total
 */
void SemiHomopolymerAlignmentPool::open(string filename, int total){
    cout<<"Load data from "<<filename<<endl;
    // clear everything
    vector<SemiHomopolymerAlignment>().swap(this->align_pool);

    // first, load in nucleotide space
    NucleotideAlignmentPool nucl_space;
    nucl_space.open(filename);
    // second, covert to homopolymer space
    unsigned rng_seed=chrono::system_clock::now().time_since_epoch().count();
    default_random_engine rng_engine(rng_seed);
    uniform_int_distribution<int> uniform_dist(1, nucl_space.align_pool.size()-total);
    int ri=uniform_dist(rng_engine)-1;

    cout<<"Use "<<total<<" out of "<<nucl_space.align_pool.size()<<" alignments"<<endl;

    vector<int> permutate(nucl_space.align_pool.size());
    for (int i=0; i<(int)nucl_space.align_pool.size(); i++){
        permutate[i]=i;
    }
    random_shuffle(permutate.begin(), permutate.end());

    for (int pi=ri; pi<ri+total; pi++){
        int i=permutate[pi];
        int other=0;
        for (int j=0; j<nucl_space.align_pool[i].raw_target.length(); j++){
            if (nucl_space.align_pool[i].raw_target[j]!='A' &&
                    nucl_space.align_pool[i].raw_target[j]!='C' &&
                    nucl_space.align_pool[i].raw_target[j]!='G' &&
                    nucl_space.align_pool[i].raw_target[j]!='T')
                other++;
        }

        if (other>0) continue;

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

void SemiHomopolymerAlignmentPool::print(string filename){
    // clear content of file
    ofstream output;
    output.open(filename, ios::out | ios::trunc);
    output.close();

    for (int i=0; i<(int)this->align_pool.size(); i++){
        align_pool[i].print(filename);
    }
}

/**
 * @brief SemiHomopolymerAlignmentPool::statistics
 * @param cycles
 */
//void SemiHomopolymerAlignmentPool::statistics(int cycles){
//   vector<SemiHomopolymerAlignmentStat> stat;
//   // do the statistics
//   statistics(stat,cycles);
//   // print out the statistics
//   for (int c=0; c<cycles; c++){
//       cout<<"[Cycle #"<<c<<"]"<<endl;
//       stat[c].print();
//   }
//}

void SemiHomopolymerAlignmentPool::statistics(int cycles){
   vector<SemiHomopolymerAlignmentStat> stat;
   // do the statistics
   statistics(stat,cycles);
   // print out the statistics
   for (int c=0; c<cycles; c++){
       SemiHomopolymerAlignmentStat s=stat[c];
       s.print(c);
   }
}

/**
 * @brief SemiHomopolymerAlignmentPool::statistics
 * @param filename
 * @param cycles
 */
void SemiHomopolymerAlignmentPool::statistics(string filename, int cycles){
    vector<SemiHomopolymerAlignmentStat> stat;
    // do the statistics
    statistics(stat,cycles);
    // print out the statistics
    ofstream ofs;
    ofs.open(filename);
    for (int c=0; c<cycles; c++){
        stat[c].print(ofs, c);
    }
    ofs.close();
}

/**
 * @brief SemiHomopolymerAlignmentPool::statistics
 * @param stat
 * @param cycles
 */
void SemiHomopolymerAlignmentPool::statistics(vector<SemiHomopolymerAlignmentStat> &stat, int cycles){
    // clear everything
    stat.clear();
    stat.resize(cycles);

    SemiHomopolymerAlignment aln1;
    vector<SemiHomopolymerAlignmentStat> aln1_stat(cycles);
    // loop over the alignments
    for (int i=0; i<(int)align_pool.size(); i++){
        // progress
        {
            char buffer[1024];
            // show progress
            int percent=(int)((i+1)*100.0/(align_pool.size()+0.));

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

            cout<< "\r"
                << "Compute statistics   "
                << "[" << bar << "] ";
            sprintf(buffer, "%3.2f", ((i+1)*100.0/(align_pool.size()+0.)));
            cout.width( 3 );
            cout<< buffer << "%     ";
            //cout<< i+1 << "/" << align_pool.size();
            cout<< flush;
        }

        // current alignment
        aln1=align_pool[i];
        // do the statistics
        aln1.statistics(aln1_stat,cycles);
        // update the overall statistician
        for (int c=0; c<cycles; c++){
            stat[c]+=aln1_stat[c];
        }
    }
    cout<<endl;
}
