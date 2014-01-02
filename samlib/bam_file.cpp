#include "bam_file.h"
#include "sam_header.h"
#include <stdio.h>
#include <iostream>
#include <map>
#include <algorithm>
#include <boost/random.hpp>
using namespace std;

BamFile::BamFile()
{
    //ctor
    bam_file_ptr=0;
}

BamFile::~BamFile()
{
    //dtor
    if (bam_file_idx) bam_index_destroy(bam_file_idx);
    if (bam_file_ptr) samclose(bam_file_ptr);
}

/** \brief To parse the BAM header field HD
 *
 * \param
 * \param
 * \return
 *
 */
void BamFile::bam_header_HD(){
    int i,n;
    const char *tags[]={"VN","SO",NULL};
    void *dict=sam_header_parse2(bam_file_ptr->header->text);
    char **tbl=sam_header2tbl_n(dict,"HD",tags,&n);

    bam_header_hd_record_n=n;
    for (i=0; i<n; i++){
        bam_header_hd_field hd;
        if (tbl[2*i]) hd.VN=tbl[2*i];
        else hd.VN="-";
        if (tbl[2*i+1]) hd.SO=tbl[2*i+1];
        else hd.SO="-";
        bam_header_hd_records.push_back(hd);
    }

    if (tbl) free(tbl);
}
/** \brief To parse the BAM header field SQ
 *
 * \param
 * \param
 * \return
 *
 */
void BamFile::bam_header_SQ(){
    int i,n;
    const char *tags[]={"SN","LN","AS","M5","SP","UR",NULL};
    void *dict=sam_header_parse2(bam_file_ptr->header->text);
    char **tbl=sam_header2tbl_n(dict,"SQ",tags,&n);

    bam_header_sq_record_n=n;
    for (i=0; i<n; i++){
        bam_header_sq_field sq;
        if (tbl[6*i])   sq.SN=tbl[6*i];   else sq.SN="-";
        if (tbl[6*i+1]) sq.LN=tbl[6*i+1]; else sq.LN="-";
        if (tbl[6*i+2]) sq.AS=tbl[6*i+2]; else sq.AS="-";
        if (tbl[6*i+3]) sq.M5=tbl[6*i+3]; else sq.M5="-";
        if (tbl[6*i+4]) sq.SP=tbl[6*i+4]; else sq.SP="-";
        if (tbl[6*i+5]) sq.UR=tbl[6*i+5]; else sq.UR="-";
        bam_header_sq_records.push_back(sq);
    }

    if (tbl) free(tbl);
}

/** \brief To parse the BAM header field RG
 *
 * \param
 * \param
 * \return
 *
 */
 void BamFile::bam_header_RG(){
    int i,n;
    const char *tags[]={"ID","CN","DS","DT","FO","KS","LB","PG","PI","PL","PU","SM",NULL};
    void *dict=sam_header_parse2(bam_file_ptr->header->text);
    char **tbl=sam_header2tbl_n(dict,"RG",tags,&n);

    bam_header_rg_record_n=n;
    for (i=0; i<n; i++){
        bam_header_rg_field rg;
        if (tbl[12*i])    rg.ID=tbl[12*i];    else rg.ID="-";
        if (tbl[12*i+1])  rg.CN=tbl[12*i+1];  else rg.CN="-";
        if (tbl[12*i+2])  rg.DS=tbl[12*i+2];  else rg.DS="-";
        if (tbl[12*i+3])  rg.DT=tbl[12*i+3];  else rg.DT="-";
        if (tbl[12*i+4])  rg.FO=tbl[12*i+4];  else rg.FO="-";
        if (tbl[12*i+5])  rg.KS=tbl[12*i+5];  else rg.KS="-";
        if (tbl[12*i+6])  rg.LB=tbl[12*i+6];  else rg.LB="-";
        if (tbl[12*i+7])  rg.PG=tbl[12*i+7];  else rg.PG="-";
        if (tbl[12*i+8])  rg.PI=tbl[12*i+8];  else rg.PI="-";
        if (tbl[12*i+9])  rg.PL=tbl[12*i+9];  else rg.PL="-";
        if (tbl[12*i+10]) rg.PU=tbl[12*i+10]; else rg.PU="-";
        if (tbl[12*i+11]) rg.SM=tbl[12*i+11]; else rg.SM="-";
        bam_header_rg_records.push_back(rg);
    }

    if (tbl) free(tbl);
 }

 /** \brief To parse the BAM header field PG
  *
  * \param
  * \param
  * \return
  *
  */
  void BamFile::bam_header_PG(){
    int i,n;
    const char *tags[]={"ID","PN","CL","PP","VN",NULL};
    void *dict=sam_header_parse2(bam_file_ptr->header->text);
    char **tbl=sam_header2tbl_n(dict,"PG",tags,&n);

    bam_header_pg_record_n=n;
    for (i=0; i<n; i++){
        bam_header_pg_field pg;
        if (tbl[5*i])   pg.ID=tbl[5*i];   else pg.ID="-";
        if (tbl[5*i+1]) pg.PN=tbl[5*i+1]; else pg.PN="-";
        if (tbl[5*i+2]) pg.CL=tbl[5*i+2]; else pg.CL="-";
        if (tbl[5*i+3]) pg.PP=tbl[5*i+3]; else pg.PP="-";
        if (tbl[5*i+4]) pg.VN=tbl[5*i+4]; else pg.VN="-";
        bam_header_pg_records.push_back(pg);
    }

    if (tbl) free(tbl);
  }

/** \brief To parse the BAM header
 *
 * \param
 * \param
 * \return
 *
 */
 void BamFile::bam_header_parser(){
    bam_header_HD();
    bam_header_SQ();
    bam_header_RG();
    bam_header_PG();
 }


 /**
  * Utilities for bam_random_retrieve
  */
 // temporary data structure for bam_random_retrieve
 typedef struct{
     int g_id, g_beg, g_end;  // the beginning position and ending position on genome
     samfile_t *bam_file_ptr;   // the pointer to the given bam file

     // items to store the read information when the region is scanned
     map<string, int> tmp_read_mark;    // mark a read as scanned and jump over next time
     vector<string> *tmp_read_pool=0;      // store the reads within the region

     // items to store the alignments
     map<string, int> *tmp_align_idx=0;
     vector<BamAlign> *tmp_align_pool=0;

     // items to store the genome content
     int tmp_genome_len;
     char *tmp_genome_seq;
     faidx_t *tmp_genome_idx;

     // filtration
     BamFilter *filter;
 }bam_random_retrieve_data_ptr;

 // temporary function for random_random_retrieve to fetch the alignments from the given bam file
static int bam_random_retrieve_fetch_func(const bam1_t *b, void *data){
    bam_plbuf_t *buf=(bam_plbuf_t*) data;
    bam_plbuf_push(b, buf);
    return 0;
}

// temporary function for bam_random_retrieve to store an read indentity into a list
static int bam_random_retrieve_read_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data){
    // temporary data information
    bam_random_retrieve_data_ptr *tmp_data=(bam_random_retrieve_data_ptr*)data;

    // loop over piled-up reads
    for (int i=0; i<n; i++){
        // current alignment
        bam_pileup1_t aln1=pl[i];
        // the leftmost and rightmost position of current alignment on genome
        int aln1_lpos=aln1.b->core.pos;
        uint32_t *aln1_cigar=bam1_cigar(aln1.b);
        int aln1_rpos=bam_calend(&aln1.b->core, aln1_cigar);
        // check the validation on the genomic position
        if (aln1_lpos<tmp_data->g_beg || aln1_rpos>tmp_data->g_end) continue;
        // check whether it pass samtools filtration
        if (aln1.b->core.flag & BAM_DEF_MASK) continue;
        // check whether it is a good mapping
        if (aln1.b->core.qual<tmp_data->filter->map_qual_thresh) continue;
        // check whether it is a long alignment
        int aln1_len=0;
        uint32_t *cigar=bam1_cigar(aln1.b);
        for (int i=0; i<aln1.b->core.n_cigar; i++){
            char op=bam_cigar_opchr(cigar[i]);
            if (op=='M' || op=='I' || op=='D') {
                aln1_len+=bam_cigar_oplen(cigar[i]);
            }
        }
        if (aln1_len<tmp_data->filter->aln_length) continue;
        // the name of current alignment
        string aln1_name=bam1_qname(aln1.b);
        if (aln1.b->core.flag & BAM_FPAIRED) {
            if (aln1.b->core.flag & BAM_FREAD1) aln1_name+=".1";
            if (aln1.b->core.flag & BAM_FREAD2) aln1_name+=".2";
        }
        // check whether the current alignment has been visited
        if (tmp_data->tmp_read_mark.count(aln1_name)==0){
            tmp_data->tmp_read_pool->push_back(aln1_name);
            tmp_data->tmp_read_mark[aln1_name]=1;
        }
    }

    return 0;
}

// temporary function for bam_random_retrieve to store an alignment into a list
static int bam_random_retrieve_align_func(uint32_t tid, uint32_t pos, int n, const bam_pileup1_t *pl, void *data){
    // temporary data information
    bam_random_retrieve_data_ptr *tmp_data=(bam_random_retrieve_data_ptr*)data;

    if (pos<tmp_data->g_beg || pos>tmp_data->g_end) return 0;

    // loop over piled-up reads
    for (int i=0; i<n; i++){
        // current alignment
        bam_pileup1_t aln1=pl[i];
        // check whether it is a good mapping
        if (aln1.b->core.qual<tmp_data->filter->map_qual_thresh) continue;
        // the name of current alignment
        string aln1_name=bam1_qname(aln1.b);
        if (aln1.b->core.flag & BAM_FPAIRED) {
            if (aln1.b->core.flag & BAM_FREAD1) aln1_name+=".1";
            if (aln1.b->core.flag & BAM_FREAD2) aln1_name+=".2";
        }
        // if current alignment is sampled, then update the record of current alignment
        if (tmp_data->tmp_align_idx->count(aln1_name)>0){
            int idx=(*tmp_data->tmp_align_idx)[aln1_name];
            // record item for current alignment
            BamAlign item=(*tmp_data->tmp_align_pool)[idx];
            // data structure for local alignment
            string local_align_read="";
            string local_align_qual="";
            string local_align_cigar="";
            string local_align_genome="";
            int local_nM=0;
            int local_nX=0;
            int local_nD=0;
            int local_nI=0;

            // data structure for the base quality scores of current alignment
            uint8_t *aln1_qual=bam1_qual(aln1.b);

            // About the definition and operation of CIGAR, refer to Samtools/bam.h for details.
            uint32_t *aln1_cigar=bam1_cigar(aln1.b);
            int aln1_cigar_n=aln1.b->core.n_cigar;

            // collect the information of the soft-clipped head
            if (aln1.is_head==1 && bam_cigar_opchr(aln1_cigar[0])=='S'){
                // TODO: pack the soft-clipped head sequence
                // pack read information
                int aln1_head_len=bam_cigar_oplen(aln1_cigar[0]);
                for (int j=aln1.qpos-aln1_head_len; j<aln1.qpos; j++){
                    local_align_read+=bam_nt16_rev_table[bam1_seqi(bam1_seq(aln1.b), j)];
                    local_align_qual+=(char)(aln1_qual[j]+33);
                }
                // pack genome information
                int tmp_genome_len;
                char *tmp_genome_data=faidx_fetch_seq(tmp_data->tmp_genome_idx, tmp_data->bam_file_ptr->header->target_name[tmp_data->g_id],
                                                      pos-aln1_head_len, pos-1, &tmp_genome_len);
                for (int j=0; j<aln1_head_len; j++){
                    if (pos+j<aln1_head_len) {
                        local_align_genome+='N';
                    }
                    else {
                        local_align_genome+=tmp_genome_data[j];
                    }
                }
                // pack alignment information
                for (int j=0; j<aln1_head_len; j++){
                    local_align_cigar+="S";
                }
            }

            // collect the information of current local alignment
            if (aln1.indel>=0){ // here is not deletion
                // HERE is match or mismatch
                // package read content
                local_align_read+=bam_nt16_rev_table[bam1_seqi(bam1_seq(aln1.b),aln1.qpos)];
                local_align_qual+=(char)(aln1_qual[aln1.qpos]+33);
                // package genome content
                local_align_genome+=tmp_data->tmp_genome_seq[pos-tmp_data->g_beg];
                // package alignment content
                if (tmp_data->tmp_genome_seq[pos-tmp_data->g_beg] == bam_nt16_rev_table[bam1_seqi(bam1_seq(aln1.b),aln1.qpos)]){
                    local_align_cigar+="M";

                    local_nM+=1;
                }else{
                    local_align_cigar+="X";

                    local_nX+=1;
                }

                // HERE has insertion
                if (aln1.indel>0){
                    for (int j=1; j<aln1.indel; j++){
                        // package read content
                        local_align_read+=bam_nt16_rev_table[bam1_seqi(bam1_seq(aln1.b),aln1.qpos+j)];
                        local_align_qual+=(char)(aln1_qual[aln1.qpos+j]+33);
                        // package genome content
                        local_align_genome+="-";
                        // package alignment content
                        local_align_cigar+="I";

                        local_nI+=1;
                    }
                }
            }else{  // here is deletion
                // package read content
                local_align_read+="-";
                local_align_qual+="!";   // "!" is 33 in ascii-table, refer to http://en.wikipedia.org/wiki/ASCII
                // package genome content
                local_align_genome+=tmp_data->tmp_genome_seq[pos-tmp_data->g_beg];
                // package alignment content
                local_align_cigar+="D";

                local_nD+=1;
            }

            // collect the information of the soft-clipped tail
            if (aln1.is_tail==1 && bam_cigar_opchr(aln1_cigar[aln1_cigar_n-1])=='S'){
                // TODO: pack the soft-clipped tail sequence
                // pack read information
                int aln1_tail_len=bam_cigar_oplen(aln1_cigar[aln1_cigar_n-1]);
                for (int j=0; j<aln1_tail_len; j++){
                    local_align_read+=bam_nt16_rev_table[bam1_seqi(bam1_seq(aln1.b), aln1.qpos+j+1)];
                    local_align_qual+=(char)(aln1_qual[aln1.qpos+j+1]+33);
                }
                // pack genome information
                int tmp_genome_len;
                char *tmp_genome_data=faidx_fetch_seq(tmp_data->tmp_genome_idx, tmp_data->bam_file_ptr->header->target_name[tmp_data->g_id],
                                                      pos+1, pos+aln1_tail_len, &tmp_genome_len);
                for (int j=0; j<aln1_tail_len; j++){
                    if (1+j>tmp_genome_len) local_align_genome+='N';
                    else local_align_genome+=tmp_genome_data[j];
                }
                // pack alignment information
                for (int j=0; j<aln1_tail_len; j++){
                    local_align_cigar+="S";
                }
            }

            // update the global item
            item.aln_read+=local_align_read;
            item.aln_read_qual+=local_align_qual;
            item.aln_genome+=local_align_genome;
            item.aln_cigar+=local_align_cigar;

            item.nM+=local_nM;
            item.nX+=local_nX;
            item.nD+=local_nD;
            item.nI+=local_nI;
            // update the mapping quality of local alignment
            if (item.aln_qual<0) {
                item.aln_qual=aln1.b->core.qual;
                if (aln1.b->core.flag & BAM_FREVERSE){
                    item.aln_strand="-";
                }else{
                    item.aln_strand="+";
                }
            }

            (*tmp_data->tmp_align_pool)[idx]=item;
        }
    }

    return 0;
}

// temporary function to generate random number
int bam_random_retrieve_random_generator(int i){
    boost::mt19937_64 no_state;
    boost::uniform_int<> rng(0, i-1);
    return rng(no_state);
}

 /**
 * @brief BamFile::bam_random_retrieve
 *             To randomly retrieve the alignments up to number in the specified region
 * @param region
 * @param number
 */
void BamFile::bam_random_retrieve(string region, int number, BamFilter filter){
    // TODO: add control to handle errors, e.g. user does not give the genome file, etc.

    vector<string> read_pool;		// all reads within region
    vector<BamAlign> align_pool;	// the alignment information of sampled reads
    map<string, int> align_idx;         // the index of an alignment in the list

    // temporary data item for the callback functions
    bam_random_retrieve_data_ptr tmp_data;
    tmp_data.bam_file_ptr=this->bam_file_ptr;
    tmp_data.tmp_read_pool=&read_pool;
    tmp_data.filter=&filter;

    // temporary buffer for bam fetch and pileup
    bam_plbuf_t *tmp_buf;

    // parse the region into items: g_id (Genome ID), g_beg (The Begion Position on Genome), g_end (The Ending Position on Genome)
    // NOTE: the user-given region starts from 1, while the Samtools API define the region as the beginning is 0.
    bam_parse_region(this->bam_file_ptr->header, region.c_str(), &tmp_data.g_id, &tmp_data.g_beg, &tmp_data.g_end);

    // #1 collect reads in the region
    tmp_buf=bam_plbuf_init(bam_random_retrieve_read_func, &tmp_data);
    bam_fetch(this->bam_file_ptr->x.bam, this->bam_file_idx, tmp_data.g_id, tmp_data.g_beg, tmp_data.g_end, tmp_buf, bam_random_retrieve_fetch_func);
    bam_plbuf_push(0, tmp_buf);

    // #2 randoly generate sample from read_pool
    random_shuffle(read_pool.begin(), read_pool.end(), bam_random_retrieve_random_generator);
    int sample_size=(number<read_pool.size())?number:read_pool.size();
    for (int i=0; i<sample_size; i++){
        BamAlign aln1;
        aln1.aln_name=read_pool[i];
        align_pool.push_back(aln1);
        align_idx[aln1.aln_name]=i;
    }

    // #3 output sampled alignments
    tmp_data.tmp_align_idx=&align_idx;
    tmp_data.tmp_align_pool=&align_pool;
    // load genome sequencing
    tmp_data.tmp_genome_idx=genome_idx;
    tmp_data.tmp_genome_seq=faidx_fetch_seq(genome_idx, bam_file_ptr->header->target_name[tmp_data.g_id], tmp_data.g_beg, tmp_data.g_end, &tmp_data.tmp_genome_len);
    // re-initialize the buffer and re-do pileup
    tmp_buf=bam_plbuf_init(bam_random_retrieve_align_func, &tmp_data);
    bam_fetch(this->bam_file_ptr->x.bam, this->bam_file_idx, tmp_data.g_id, tmp_data.g_beg, tmp_data.g_end, tmp_buf, bam_random_retrieve_fetch_func);
    bam_plbuf_push(0, tmp_buf);
    // output the randomly sampled alignments
    for (int i=0; i<align_pool.size(); i++){
        BamAlign item=align_pool[i];
        item.print();
    }

    // #4 destroy the pileup buffer
    bam_plbuf_destroy(tmp_buf);
}

/** \brief To print out the generic information of BAM file
 *
 * \param
 * \param
 * \return
 *
 */
void BamFile::bam_generic_info(){
    path bam_file_path(bam_file_name);
    cout<<"Directory: "<<bam_file_path.parent_path()<<endl;
    cout<<"Filename:  "<<bam_file_path.filename()<<endl;
    cout<<"[Header Field HD]"<<endl;
    for (int i=0; i<bam_header_hd_record_n; i++){
        cout<<"VN:"<<bam_header_hd_records[i].VN<<"  SO:"<<bam_header_hd_records[i].SO<<endl;
    }
    cout<<"[Header Field SQ]"<<endl;
    for (int i=0; i<bam_header_sq_record_n; i++){
        cout<<"SN:"<<bam_header_sq_records[i].SN<<"  LN:"<<bam_header_sq_records[i].LN<<"  "
            <<"AS:"<<bam_header_sq_records[i].AS<<"  M5:"<<bam_header_sq_records[i].M5<<"  "
            <<"SP:"<<bam_header_sq_records[i].SP<<"  UR:"<<bam_header_sq_records[i].UR<<endl;
    }
    cout<<"[Header Field RG]"<<endl;
    for (int i=0; i<bam_header_rg_record_n; i++){
        cout<<"ID:"<<bam_header_rg_records[i].ID<<"  CN:"<<bam_header_rg_records[i].CN<<"  "
            <<"DS:"<<bam_header_rg_records[i].DS<<"  DT:"<<bam_header_rg_records[i].DT<<"  "
            <<"FO:"<<bam_header_rg_records[i].FO<<"  KS:"<<bam_header_rg_records[i].KS<<"  "
            <<"LB:"<<bam_header_rg_records[i].LB<<"  PG:"<<bam_header_rg_records[i].PG<<"  "
            <<"PI:"<<bam_header_rg_records[i].PI<<"  PL:"<<bam_header_rg_records[i].PL<<"  "
            <<"PU:"<<bam_header_rg_records[i].PU<<"  SM:"<<bam_header_rg_records[i].SM<<endl;
    }
    cout<<"[Header Field PG]"<<endl;
    for (int i=0; i<bam_header_pg_record_n; i++){
        cout<<"ID:"<<bam_header_pg_records[i].ID<<"  PN:"<<bam_header_pg_records[i].PN<<"  "
            <<"CL:"<<bam_header_pg_records[i].CL<<"  PP:"<<bam_header_pg_records[i].PP<<"  "
            <<"PP:"<<bam_header_pg_records[i].PP<<"  VN:"<<bam_header_pg_records[i].VN<<endl;
    }
}
