#ifndef SEMI_HOMO_GHMM_ORDER1_H
#define SEMI_HOMO_GHMM_ORDER1_H

#include "semi_homo_align.h"

namespace SemiHomopolymerAlignmentSpace{
    enum AlignmentState{A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,Ag15,    // g15 is "greater than 15"
                        C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,Cg15,    // g15 is "greater than 15"
                        G1,G2,G3,G4,G5,G6,G7,G8,G9,G10,G11,G12,G13,G14,G15,Gg15,    // g15 is "greater than 15"
                        T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14,T15,Tg15,    // g15 is "greater than 15"
                        I,B,E};
    const int ALIGNMENTSTATSIZE=4*16+1+2;
}

class SemiHomopolymerGHMMOrder1
{
public:
    SemiHomopolymerGHMMOrder1();
};

#endif // SEMI_HOMO_GHMM_ORDER1_H
