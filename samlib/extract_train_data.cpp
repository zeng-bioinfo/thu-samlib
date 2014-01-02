#include "OptionPrinter.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <stdio.h>
#include "sam.h"
#include "bam_file.h"

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace


int main(int argc, char** argv)
{
  try
  {
    std::string appName = boost::filesystem::basename(argv[0]);
    std::string bam_filename, genome_filename;
    std::string region;
    int number=1000;
    int block_size=1000;
    int map_qual_thresh=10;
    int aln_length=100;


    // Define and parse the program options
    namespace po = boost::program_options;
    po::options_description desc("Options");
    desc.add_options()
        ("bamfile,b", po::value<std::string>(&bam_filename)->required(), "The filename of BAM file")
        ("genome,g",po::value<std::string>(&genome_filename)->required(), "The filename of genome sequence")
        ("region,r",po::value<std::string>(&region)->required(), "Region specification, e.g. chr20:1000000-1500000")
        ("number,n",po::value<int>(&number), "The number of alignments would be extracted out from the specified region [1000]")
        ("block_size,z",po::value<int>(&block_size), "It will divide the region into segments of block_sizes, then extract from each segment evenly [1000]")
        ("map_qual_thres,m",po::value<int>(&map_qual_thresh), "The threshod of the mapping quality score [10]")
        ("aln_length,l",po::value<int>(&aln_length), "The threshold of the alignment length [100]")
        ("help,h", "Print help messages");

    po::variables_map vm;

    try
    {
      po::store(po::command_line_parser(argc, argv).options(desc).run(),vm); // throws on error

      // --help option
      if ( vm.count("help")  )
      {
        std::cout << "This program is to extract the alignments from the BAM file"
                  << " to be the traning data"
                  << std::endl << std::endl;
        rad::OptionPrinter::printStandardAppDesc(appName,
                                                 std::cout,
                                                 desc);
        return SUCCESS;
      }

      po::notify(vm); // throws on error, so do after help in case
                      // there are any problems
    }
    catch(boost::program_options::required_option& e)
    {
      rad::OptionPrinter::formatRequiredOptionError(e);
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      rad::OptionPrinter::printStandardAppDesc(appName,
                                               std::cout,
                                               desc);
      return ERROR_IN_COMMAND_LINE;
    }
    catch(boost::program_options::error& e)
    {
      std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
      rad::OptionPrinter::printStandardAppDesc(appName,
                                               std::cout,
                                               desc);
      return ERROR_IN_COMMAND_LINE;
    }

    // using Samtools API to parse the region specification
    int g_id, g_beg, g_end;
    std::string g_name;

    samfile_t *bam_ptr=samopen(bam_filename.c_str(), "rb", NULL);
    bam_parse_region(bam_ptr->header, region.c_str(), &g_id, &g_beg, &g_end);
    g_name=bam_ptr->header->target_name[g_id];
    samclose(bam_ptr);

    // division
    int block_num=(g_end-g_beg)/block_size;
    int block_sample_size=number/block_num;

    // sampling
    int i;
    char buffer[1024];
    BamFilter filter;
    filter.map_qual_thresh=map_qual_thresh;
    filter.aln_length=aln_length;
    BamFile bamfile(bam_filename, genome_filename);
    for (i=0; i<block_num-1; i++){
        sprintf(buffer,"%s:%d-%d",g_name.c_str(),g_beg+i*block_size,g_beg+i*block_size+block_size-1);
        bamfile.bam_random_retrieve(std::string(buffer),block_sample_size,filter);
    }
    sprintf(buffer,"%s:%d-%d",g_name.c_str(),g_beg+i*block_size,g_beg+i*block_size+block_size-1);
    bamfile.bam_random_retrieve(std::string(buffer),block_sample_size+number%block_num,filter);
//*/


  }
  catch(std::exception& e)
  {
    std::cerr << "Unhandled Exception reached the top of main: "
              << e.what() << ", application will now exit" << std::endl;
    return ERROR_UNHANDLED_EXCEPTION;

  }

  return SUCCESS;

} // main
