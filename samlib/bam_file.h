#ifndef BAMFILE_H
#define BAMFIEL_H
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include "sam.h"
using namespace std;
using namespace boost::filesystem;

/**
 * Data structure for a BAM alignment.
 * This class is for single-end read.
 */
class BamAlign{
	string aln_name;			// read id
	string aln_read;			// read symbols
	vector<int> aln_read_qual;	// read quality scores
	string aln_genome;			// reference symbols
	string aln_cigar;			// alignment types
	int aln_qual;				// alignment quality score

	BamAlign& operator=(const BamAlign& _aln){
		this->aln_name=_aln.aln_name;
		this->aln_read=_aln.aln_read;
		for (int i=0; i<_aln.aln_read_qual.size(); i++)
			this->aln_read_qual.push_back(_aln.aln_read_qual[i]);
		this->aln_genome=_aln.aln_genome;
		this->aln_cigar=_aln.aln_cigar;
		this->aln_qual=_aln.aln_qual;
		return *this;
	}
	/**
	 * Default constructor
	 */
	BamAlign(){ ; }
	/**
	 * Copy constructor
	 * @param _aln
	 */
	BamAlign(const BamAlign& _aln){
		this->aln_name=_aln.aln_name;
		this->aln_read=_aln.aln_read;
		for (int i=0; i<_aln.aln_read_qual.size(); i++)
			this->aln_read_qual.push_back(_aln.aln_read_qual[i]);
		this->aln_genome=_aln.aln_genome;
		this->aln_cigar=_aln.aln_cigar;
		this->aln_qual=_aln.aln_qual;
	}
	/**
	 * Default deconstructor
	 */
	~BamAlign(){ ; }
};

/**
 * Data structure for a BAM file
 */
class BamFile
{




//--------------------------------------------------------------------------------------------------------------------

    /**< data structures and operations for BAM header */
    public:
        /**< data structure for BAM header field HD */
        struct bam_header_hd_field{
            string VN;  // format version
            string SO;  // sorting order of alignments. Valid values: unknown, unsorted, queryname and coordinate
        };
        /**< data structure for BAM header field SQ */
        struct bam_header_sq_field{
            string SN;  // sequence name
            string LN;  // sequence length
            string AS;  // assembly identifier
            string M5;  // MD5 checksum
            string SP;  // species
            string UR;  // URI of the sequence
        };
        /**< data structure for BAM header field RG */
        struct bam_header_rg_field{
            string ID;  // read group identifier
            string CN;  // name of sequencing center producing the read
            string DS;  // description
            string DT;  // date the run was produced
            string FO;  // flow order
            string KS;  // the array of nucleotide bases that correspond to the key sequence of each read
            string LB;  // library
            string PG;  // program used for processing the read group
            string PI;  // predicted median insert size
            string PL;  // platform, valid values: CAPILLARY,LS454,ILLUMINA,SOLID,HELICOS,IONTORRENT and PACBIO
            string PU;  // platform unit (e.g. lane for Illumina or slide for SOLiD)
            string SM;  // sample
        };
        /**< data structure for BAM header field PG */
        struct bam_header_pg_field{
            string ID;  // program record identifier
            string PN;  // program name
            string CL;  // command line
            string PP;  // previous @PG-ID
            string VN;  // program version
        };
    public:
        /**< record for BAM header field HD */
        int bam_header_hd_record_n;
        vector<bam_header_hd_field> bam_header_hd_records;
        /**< record for BAM header field SQ */
        int bam_header_sq_record_n;
        vector<bam_header_sq_field> bam_header_sq_records;
        /**< record for BAM header field RG */
        int bam_header_rg_record_n;
        vector<bam_header_rg_field> bam_header_rg_records;
        /**< record for BAM header field PG */
        int bam_header_pg_record_n;
        vector<bam_header_pg_field> bam_header_pg_records;
    private:
        /**< top-level routine to parse all BAM header fields */
        void bam_header_parser();
        /**< routine to parse BAM header field HD */
        void bam_header_HD();
        /**< routine to parse BAM header field SQ */
        void bam_header_SQ();
        /**< routine to parse BAM header field RG */
        void bam_header_RG();
        /**< routine to parse BAM header field PG */
        void bam_header_PG();




//--------------------------------------------------------------------------------------------------------

    /**< data structures for BAM object */
    public:
        /**< the name of BAM file */
        string bam_file_name;
        /**< the pointer to BAM file */
        samfile_t *bam_file_ptr;
    public:
        /**< initialize the pointer */
        void init() {
            path bam_file_path(bam_file_name);
            if (!exists(bam_file_path)) {
                cerr<<bam_file_path<<" does not exist\n";
                exit(EXIT_FAILURE);
            }
            bam_file_ptr=samopen(bam_file_name.c_str(),"rb",NULL);
            bam_header_parser();
        }
        void init(string fn){
            bam_file_name=fn;
            init();
        }



//-------------------------------------------------------------------------------------------------
        /**
         * To randomly retrieve alignments of 'number' in the specified 'region'.
         * @param region
         * @param number
         */
         void bam_random_retrieve(string region, int number);


//-------------------------------------------------------------------------------------------------

    /**< verbose mode */
    public:
        /**< print out the generic information of BAM file */
        void bam_generic_info();




//--------------------------------------------------------------------------------------------------

    /**< constructor and destructor */
    public:
        BamFile();
        BamFile(string fn):bam_file_name(fn){
            init();
        }
        virtual ~BamFile();
    protected:
    private:
};


#endif // MYBAM_H
