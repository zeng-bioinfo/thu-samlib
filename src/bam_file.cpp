#include "MyBam.h"
#include <samtools/sam_header.h>
#include <iostream>
using namespace std;

BamFile::BamFile()
{
    //ctor
    bam_file_ptr=0;
}

BamFile::~BamFile()
{
    //dtor
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
