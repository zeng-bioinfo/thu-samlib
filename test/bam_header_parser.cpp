#include <string>
#include <iostream>
#include <samtools/sam.h>
#include <samtools/sam_header.h>
#include "MyBam.h"

using namespace std;

int main(int argc, char** argv)
{

    string fn=argv[1];
    BamFile bamfile(fn);
    bamfile.bam_generic_info();

    return 0;
}
