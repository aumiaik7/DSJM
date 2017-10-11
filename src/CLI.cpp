#include "CLI.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

CLI::ordering_method CLI::parseMethod(char *string)
{
    if(string != NULL && strlen(string) > 0 )
    {
        if(strcmp("lfo",string) == 0)
            return LFO;
        else if ( strcmp ("ido" , string) == 0)
            return IDO;
        else if (strcmp ( "slo" , string) == 0)
            return SLO;
        else if (strcmp ( "rlf" , string) == 0)
            return RLF;
        else if (strcmp ( "sdo" , string) == 0)
            return SDO;
        else if (strcmp ( "profit" , string) == 0)
            return PROFIT;
        else if (strcmp ( "profit_rlf", string) == 0)
            return PROFIT_RLF;
        else if (strcmp ( "rlf_slo" , string) == 0 )
            return RLF_SLO;
        else if (strcmp ( "rlf_ido", string) == 0)
            return RLF_IDO;
        else if (strcmp ( "rlf_lfo", string) == 0)
            return RLF_LFO;
        else if (strcmp ( "slo_rlf", string) == 0)
            return SLO_RLF;
        else if (strcmp ( "ido_rlf", string) == 0)
            return IDO_RLF;
        else if (strcmp ( "lfo_rlf", string) == 0)
            return LFO_RLF;
  	else if (strcmp ( "exact", string) == 0)
            return EXACT;
        else
            return UNKNOWN;
    }
    else
        return UNKNOWN;
}

CLI::approximation_method CLI::parseApproximationMethod(char *string)
{
    if(string != NULL && strlen(string) > 0 )
    {
        if(strcmp("fd",string) == 0)
            return FD;
        else if ( strcmp ("cfd" , string) == 0)
            return CFD;
        else if (strcmp ( "dd" , string) == 0)
            return DD;
        else
            return a_UNKNOWN;
    }
    else
        return a_UNKNOWN;
}

string CLI::approximation_method_to_string(CLI::approximation_method o)
{
    if(o == FD)
        return "FD";
    else if( o == CFD)
        return "CFD" ;
    else if (o == DD)
        return "DD";
    else
        return "UNKNOWN";
}

string CLI::method_to_string(ordering_method o)
{
    if(o == LFO)
        return "LFO";
    else if( o == IDO)
        return "IDO" ;
    else if (o == RLF)
        return "RLF";
    else if (o == SLO)
        return "SLO";
    else if (o == SDO)
        return "SDO";
    else if (o == PROFIT)
        return "PROFIT";
    else if (o == PROFIT_RLF)
        return "PROFIT_RLF";
    else if ( o == RLF_SLO )
        return "RLF_SLO";
    else if (o == RLF_IDO)
        return "RLF_IDO";
    else if ( o == RLF_LFO)
        return "RLF_LFO";
    else if ( o == SLO_RLF)
        return "SLO_RLF";
    else if ( o == IDO_RLF)
        return "IDO_RLF";
    else if ( o == LFO_RLF)
        return "LFO_RLF";
     else if ( o == EXACT)
        return "EXACT";	
    else
        return "UNKNOWN";
}
