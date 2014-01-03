#include "nucl_align.h"
#include "OptionPrinter.h"
#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>
#include <iostream>
#include <string>

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
      int cycles=1;


      // Define and parse the program options
      namespace po = boost::program_options;
      po::options_description desc("Options");
      desc.add_options()
          ("filename,f", po::value<std::string>(&filename)->required(), "The filename of the alignments")
          ("mode,m", po::value<std::string>(&mode), "Valid option values: nucleotide or homopolymer [nucleotide]")
          ("cycles,c",po::value<int>(&cycles), "The number of virtual sequencing cycles [1]")
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

        // --mode option
        if (vm.count("mode")){
            if (mode!="nucleotide" || mode!="homopolymer"){
                throw po::validation_error(po::validation_error::invalid_option_value, "mode");
            }
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
      catch (boost::program_options::validation_error& e){
          std::cerr << "ERROR: " << e.what() << std::endl << std::endl;
          rad::OptionPrinter::printStandardAppDesc(appName,
                                                   std::cout,
                                                   desc);
          return ERROR_IN_COMMAND_LINE;
      }

      // do the thing
      if (mode=="nucleotide"){  // nucleotide space
        NucleotideAlignmentPool align_pool;
        align_pool.open(filename);
      }else{

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
