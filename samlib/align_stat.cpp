#include "nucl_align.h"
#include "semi_homo_align.h"
#include "OptionPrinter.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>
#include <vector>

namespace
{
  const size_t ERROR_IN_COMMAND_LINE = 1;
  const size_t SUCCESS = 0;
  const size_t ERROR_UNHANDLED_EXCEPTION = 2;

} // namespace

int main(int argc, char **argv){
    try
    {
      std::string appName = boost::filesystem::basename(argv[0]);
      std::string filename,mode="nucleotide";
      std::string outfile;
      int cycles=1;
      bool plt_dot_matrix=false;


      // Define and parse the program options
      namespace po = boost::program_options;
      po::options_description desc("Options");
      desc.add_options()
          ("filename,f", po::value<std::string>(&filename)->required(), "The filename of the alignments")
          ("output,o", po::value<std::string>(&outfile), "The filename of the output")
          ("mode,m", po::value<std::string>(&mode), "Valid option values: nucleotide or homopolymer [nucleotide]")
          ("cycles,c",po::value<int>(&cycles), "The number of virtual sequencing cycles [1]")
          ("dot_matrix,d","Toggle on the dot-matrix output")
          ("help,h", "Print help messages");

      po::variables_map vm;

      try
      {
        po::store(po::command_line_parser(argc, argv).options(desc).run(),vm); // throws on error

        // --help option
        if ( vm.count("help")  )
        {
          std::cout << "This program is to calculate the statistics from the alignments"
                    << " in the nucleotide space or homopolymer space"
                    << std::endl << std::endl;
          rad::OptionPrinter::printStandardAppDesc(appName,
                                                   std::cout,
                                                   desc);
          return SUCCESS;
        }


        po::notify(vm); // throws on error, so do after help in case
                        // there are any problems

        // --mode option
        if (vm.count("mode")){
            if (mode!="nucleotide" && mode!="homopolymer"){
                throw po::validation_error(po::validation_error::invalid_option_value, string("mode:"+mode).c_str());
            }
        }
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

      if (vm.count("dot_matrix")){
          plt_dot_matrix=true;
      }

      // do the thing
      if (!plt_dot_matrix){
          if (mode=="nucleotide"){  // nucleotide space
              NucleotideAlignmentPool align_pool;
              align_pool.open(filename);
              if(vm.count("output")){
                  align_pool.statistics(outfile, cycles);
              }else{
                  align_pool.statistics(cycles);
              }
          }else{                    // homopolymer space
              SemiHomopolymerAlignmentPool align_pool;
              align_pool.open(filename);
              align_pool.statistics(outfile, cycles);
          }
      }else{
          NucleotideAlignmentTool dot_matrix_proxy;
          NucleotideAlignmentPool align_pool;
          align_pool.open(filename);
          for (int i=0; i<(align_pool.align_pool.size()); i++){
              dot_matrix_proxy.print_dot_matrix(align_pool.align_pool[i]);
          }
      }


    }
    catch(std::exception& e)
    {
      std::cerr << "Unhandled Exception reached the top of main: "
                << e.what() << ", application will now exit" << std::endl;
      return ERROR_UNHANDLED_EXCEPTION;

    }

    return SUCCESS;
}
