#include <iostream>
#include <string>
#include "bam_file.h"
using namespace std;

int main(int argc, char** argv){
    // bam filename
    string bam_filename=argv[1];
    string genome_filename=argv[2];

    // BamFile object
    BamFile bamfile(bam_filename, genome_filename);
    // genomic region
    string region="chr20:20000000-20000500";
    // filtration
    BamFilter filter;
    // retrieve alignments within the region
    bamfile.bam_random_retrieve(region, 10, filter);

    return 0;
}
