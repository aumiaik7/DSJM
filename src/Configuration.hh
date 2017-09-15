#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include "CLI.h"

#include <string>


class Configuration
{
public:
    Configuration();
    virtual ~Configuration();

    CLI::ordering_method oMethod;
    char* inputfile; /* Matrix Description in Matrix Market format */
    bool statistics_only;
    bool is_cseg;
    char* partitionFile;
    bool extend;
    bool verify;
    bool inputfile_given;
    bool io_at_end;

    string exeName;
    string commandLine;

    bool is_symmetric;
    bool is_pattern;

    bool is_load_each_row_as_partition;

    int M;
    int N;
    int nz;
    int NZ;
    int tieBrkDsat;
};

#endif
