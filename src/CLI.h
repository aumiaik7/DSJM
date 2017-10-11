#ifndef CLI_H
#define CLI_H
#include <string>
using namespace std;
namespace CLI
{
    enum ordering_method {
        IDO,
        IDO_RLF,
        LFO,
        LFO_RLF,
        NNZ_RLF,
        PROFIT,
        PROFIT_RLF,
        RLF,
        RLF_IDO,
        RLF_LFO,
        RLF_SLO,
        SDO,
        SLO,
        SLO_RLF,
	    EXACT,
        UNKNOWN,

    };

    enum approximation_method{
        FD, // Forward Differencing
        CFD, // Compressed Forward Differencing
        DD, // Direct Differentiation
        a_UNKNOWN, // Unknown
    };
    ordering_method parseMethod(char *string);
    approximation_method parseApproximationMethod(char *string);
    string method_to_string(ordering_method o);
    string approximation_method_to_string(approximation_method o);
}

#endif
