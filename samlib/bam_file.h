#ifndef BAMFILE_H
#define BAMFIEL_H
#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include "sam.h"
#include "faidx.h"
using namespace std;
using namespace boost::filesystem;

/**
 * Data structure for a BAM alignment.
 * This class is for single-end read.
 */
class BamAlign{
	public:
		string aln_name;			// read id
		string aln_read;			// read symbols
        string aln_read_qual;	    // read quality scores
		string aln_genome;			// reference symbols
		string aln_cigar;			// alignment types
		int aln_qual;				// alignment quality score
        string aln_strand;          // alignment direction

	public:
		BamAlign& operator=(const BamAlign& _aln){
			this->aln_name=_aln.aln_name;
			this->aln_read=_aln.aln_read;
			this->aln_read_qual=_aln.aln_read_qual;
			this->aln_genome=_aln.aln_genome;
			this->aln_cigar=_aln.aln_cigar;
			this->aln_qual=_aln.aln_qual;
            this->aln_strand=_aln.aln_strand;
			return *this;
		}
		void print(){
            cout<<"> "<<this->aln_name<<" mq:"<<this->aln_qual<<" strand:"<<this->aln_strand<<endl
				<<this->aln_read<<endl
				<<this->aln_read_qual<<endl
				<<this->aln_genome<<endl
				<<this->aln_cigar<<endl;
		}
	public:
		/**
		 * Default constructor
		 */
        BamAlign(){
            this->aln_name="";
            this->aln_genome="";
            this->aln_read="";
            this->aln_read_qual="";
            this->aln_cigar="";
            this->aln_qual=-1;
            this->aln_strand="";
        }
		/**
		 * Copy constructor
		 * @param _aln
		 */
		BamAlign(const BamAlign& _aln){
			this->aln_name=_aln.aln_name;
			this->aln_read=_aln.aln_read;
			this->aln_read_qual=_aln.aln_read_qual;
			this->aln_genome=_aln.aln_genome;
			this->aln_cigar=_aln.aln_cigar;
			this->aln_qual=_aln.aln_qual;
            this->aln_strand=_aln.aln_strand;
		}
		/**
		 * Default deconstructor
		 */
		virtual ~BamAlign(){ ; }
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
        /**
         * @brief bam_file_name
         */
        string bam_file_name;
        /**
         * @brief bam_file_ptr
         */
        samfile_t *bam_file_ptr;
        /**
         * @brief bam_file_idx
         */
        bam_index_t *bam_file_idx;
    private:
        /**
         * @brief init    initialize the pointer
         */
        void init() {
            // set up the content about BAM file
            path bam_file_path(bam_file_name);
            if (!exists(bam_file_path)) {
                cerr<<bam_file_path<<" does not exist\n";
                exit(EXIT_FAILURE);
            }
            bam_file_ptr=samopen(bam_file_name.c_str(),"rb",NULL);
            bam_file_idx=bam_index_load(bam_file_name.c_str());
            bam_header_parser();
            // set up the content about GENOME file
            if (!genome_file_name.empty()) genome_idx=fai_load(genome_file_name.c_str());
        }
    public:
        void init(string bn, string gn){
            bam_file_name=bn;
            genome_file_name=gn;
            init();
        }



    public:
        /**
         * @brief genome_file_name
         */
        string genome_file_name;
        /**
         * @brief genome_idx
         */
        faidx_t *genome_idx;


//-------------------------------------------------------------------------------------------------
    public:
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
        BamFile(string bn, string gn):bam_file_name(bn), genome_file_name(gn){
            init();
        }
        virtual ~BamFile();
    protected:
    private:
};


#endif // MYBAM_H
