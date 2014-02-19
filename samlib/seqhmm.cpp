#include "semi_homo_align.h"
#include "semi_homo_ghmm_order1.h"
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
        std::string traindata;
        std::string output;
        std::string output_binary;
        std::string mode="re-alignment";
        std::string model_config;
        int markov_order=1;
        int cycles=1;
        int band=10;
        int number=1000;


        // Define and parse the program options
        namespace po = boost::program_options;
        po::options_description desc("Options");
        desc.add_options()
                ("mode,m", po::value<std::string>(&mode), "Running mode: train or re-alignment [re-alignment]")
                ("number", po::value<int>(&number), "The number of alignments used to train model [1000]")
                ("Markov-chain-order", po::value<int>(&markov_order), "The order of the hidden Markov chain [1]")
                ("cycles,c", po::value<int>(&cycles), "The number of virtual sequencing cycles [1]")
                ("band", po::value<int>(&band), "Specify alignment width by 'band' percent [10]")
                ("data", po::value<std::string>(&traindata), "The filename of original alignment data")
                ("model-config", po::value<std::string>(&model_config), "The filename of model configuration")
                ("output,o", po::value<std::string>(&output), "The filename of output data (for debug)")
                ("output-binary,b", po::value<std::string>(&output_binary), "The filename of output binary data")
                ("help,h", "Print help messages");

        po::variables_map vm;

        try
        {
            po::store(po::command_line_parser(argc, argv).options(desc).run(),vm); // throws on error

            // --help option
            if ( vm.count("help") || vm.empty() )
            {
                std::cout << "The generalized hidden Markov model to model the statistical characteristics of"
                          << "next-generation sequencing technologies, such as Illumina, Ion Torrent, 454, and"
                          << "Pacific Biosciences\n\n"
                          << "This program is to train the parameters of the generalized hidden Markov model"
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
                if (mode!="train" && mode!="re-alignment"){
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


        // do the thing
        if (mode=="train"){   // train the generalized hidden Markov model
            if (markov_order==1){
                SemiHomopolymerGHMMOrder1 ghmm(cycles);

                SemiHomopolymerAlignmentPool data;
                if (vm.count("number")){
                    data.open(traindata, number);
                }else{
                    data.open(traindata);
                }

                // train the model
                ghmm.parameter_estimate(data, band, 100, 1e-5);

                if (!output.empty()){
                    ghmm.parameter_print(output);
                }

                if (!output_binary.empty()){
                    ghmm.parameter_print_binary(output_binary);
                }

                if (output.empty() && output_binary.empty()){
                    ghmm.parameter_print();
                }
            }
        }

        if (mode=="re-alignment"){
            SemiHomopolymerGHMMOrder1 ghmm;
            ghmm.load_model_configure(model_config);

            SemiHomopolymerAlignmentPool data;
            if (vm.count("number")){
                data.open(traindata, number);
            }else{
                data.open(traindata);
            }

            for (int i=0; i<(int)data.align_pool.size(); i++){
                cout<<"\r" "Compute re-alignment "<<i+1<<" out of "<<data.align_pool.size()<<flush;

                SemiHomopolymerAlignment aln1;
                ghmm.compute_banded_alignment(data.align_pool[i].align_name,
                                              data.align_pool[i].raw_target,
                                              data.align_pool[i].raw_query,
                                              data.align_pool[i].raw_quality,
                                              data.align_pool[i].query_strand,
                                              aln1, band);
                data.align_pool[i]=aln1;
            }
            cout<<endl;

            if (!output.empty()){
                data.print(output);
            }else{
                data.print();
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

