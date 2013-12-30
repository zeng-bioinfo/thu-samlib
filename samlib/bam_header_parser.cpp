#include <string>
#include <iostream>
#include "sam.h"
#include "sam_header.h"
#include "bam_file.h"

using namespace std;

int main(int argc, char** argv)
{

    string bn=argv[1];
    string gn=argv[2];
    BamFile bamfile(bn, gn);
    bamfile.bam_generic_info();

    return 0;
}
