#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <sys/times.h>
#include <sys/resource.h>
#include <cstring>
#include <fstream>
#include <string>
#include "Timer.h"

#include "CLI.h"
#include "gcolor.hh"
#include "mmio.h"
#include "Matrix.hh"
#include "Result.hh"
#include "Configuration.hh"
#include "driver.hh"

#include "Utility.h"

using namespace std;

/**
 *
 * The following optString is an important part of the Command Line Interface.
 * -m lfo -h -? -i inputfile.mtx
 */
static const char *optString = "m:i:h?sp:evrqct:";
// -c The Command Line itself, we want to print the command line
// -m rlf/ido/slo !Name of the Ordering Algorithm!
// -i matrix.mtx !Name of the Matrix File!
// -h !Help!
// -s !Statistics Only!
// -e !Extend The Matrix!
// -v !Generate Results for Verification!
// -r !IO at end only!
// -p partition.mtx
// -q load each row as partition.


/**
 * This function is called when the application is called with the display help
 * command Line options.
 */
void display_usage(void)
{
  puts("Graph Coloring Driver Program");
  puts("Usage: ");
  puts("        -m rlf/ido/sdo/slo/lfo/exact !Name of the Ordering Algorithm!");
  puts("        -i matrix.mtx !Name of the matrix File!");
  puts("        -h !Help!");
  puts("        -s !Statistics Only!");
  puts("        -e !Extend the Matrix!");
  puts("        -v !Generate Results for Verification!");
  puts("        -q !Load Each row as a partition!");
  puts("        -c !The Command Line!");
  puts("        -t !Only for DSATUR-exact. DSATUR-exact tie bracking method(1/2/3/4)!");
  /* ... */
  exit(EXIT_FAILURE);
}
/**
 * \brief
 * Description:
 * This function runs a loop to process all the command line
 * parameters and then set the corresponding members in the global
 * struct variable "configuration".
 */
Configuration* processCommandLineParameters(int argc, char **argv, Configuration* configuration )
{
    int opt  = 0;

    /**
     * Initialize the configuration before we get to work. */
    configuration->oMethod = CLI::SLO;
    configuration->statistics_only = false;
    configuration->is_cseg = false;
    configuration->extend = false;
    configuration->verify = false;
    configuration->inputfile_given = false;
    configuration->is_load_each_row_as_partition = false ;
    configuration->tieBrkDsat = 1;


    /* We are starting parsing the command Line arguments with getopt
     * method
     */

    opt = getopt(argc,argv, optString);
    while ( opt != -1 )
    {
        switch( opt )
        {
        case 'c':
            // The Command Line itself.
            // std::string s(optarg);
            configuration->commandLine = MatrixUtility::createString(optarg);
        case 'm' :
            // Select the method for ordering
            // between SLO, IDO and LFO
            // RLF is handled by a separate program.

            // But we might merge them in a single program. in near future.
            // configuration->oMethod = optarg;
            configuration->oMethod = CLI::parseMethod(optarg);
            if(configuration->oMethod == CLI::LFO)
                // cout << "You have selected LFO" << endl;
                ;
            else if( configuration->oMethod == CLI::SLO )
                // cout << "You have selected SLO" << endl;
                ;
            else if ( configuration->oMethod == CLI::IDO )
                // cout << "You have selected IDO " << endl;
                ;
            else if( configuration->oMethod == CLI::RLF  )
                // cout << "You have selected RLF" << endl;
                ;
            else if(configuration->oMethod == CLI::SDO)
                ;
	    else if (configuration->oMethod == CLI::SLO_RLF)
                ;

            else if (configuration->oMethod == CLI::RLF_SLO)
                ;
	    else if (configuration->oMethod == CLI::EXACT)
                ;	

            else
                cout << "Unknown method selected for ordering " << endl;
            break;
        case 'i':
            /**
             * \todo Should we change the struct configuration to a class
             * so that it can better serve in case dealing with memory handling ?
             * So we got something with the -i option.
             * This is a character string and most probably the name
             * of the file containing the matrix description in matrix
             * market format. We should store this in our configuration struct.
             */
            if(optarg != NULL && strlen(optarg ) > 0)
            {
                if(configuration->inputfile != NULL)
                    delete[] configuration->inputfile ;
                configuration->inputfile = new char[strlen(optarg) + 1];
                strcpy(configuration->inputfile, optarg);
                // cout << "We have copied it to the configuration - " <<
                // configuration->inputfile << endl;
                configuration->inputfile_given = true;
            }

            break;
        case 'h': /* fall-through is intentional */
        case '?':
            display_usage();
        break;
        case 's':
            configuration->statistics_only = true;
            break;

        case 'p':
            if(optarg != NULL && strlen(optarg) > 0)
            {
                delete[] configuration->partitionFile;
                configuration->partitionFile = new char[strlen(optarg) + 1] ;
                strcpy(configuration->partitionFile, optarg);

                // Set the boolean variable indicating that this is a
                // case of CSeg matrix.
                configuration->is_cseg = true;
            }
            break;
        case 'e': // Extend-Only
            configuration->extend = true;
            break;
        case 'v': // Verify
            configuration->verify = true;
            break;
        case 'r': // IO At End Only, not when clock is running
            configuration->io_at_end = false;
            break;
        case 'q':
            configuration->is_load_each_row_as_partition = true;
            configuration->is_cseg = true;
            configuration->partitionFile = NULL;
            break;
        case 't': 
            if(configuration->oMethod == CLI::EXACT && atoi(optarg) <=4)
            {
                configuration->tieBrkDsat = atoi(optarg);
            }
            else
            {
                cout<<"-t option is only for dsatur-exact tie braking scheme."<<endl;  
                display_usage();
            }
            
            break;
        default:
            break;
        }

        opt = getopt(argc,argv,optString);
    }

    return configuration;

}

int main(int argc, char *argv[])
{
  Configuration *configuration = new Configuration();
  configuration->exeName  = string(argv[0]);

  if (argc < 2 )
    {
      fprintf (stderr, "Usage: %s -i [martix-market-filename] -m [slo/ido/rlf/lfo/sdo/exact]\n", argv[0]);
      exit(1);
    }
  /**
   *
   *
   * Process Command Line Parameters
   * ----------------------------------------
   *
   */
  configuration = processCommandLineParameters(argc,argv, configuration);

  run(configuration);
  delete configuration;
  return 0;
}


