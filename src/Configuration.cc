#include "Configuration.hh"

/**
 * @Description : Constructor
 * @Creation Date: 05.08.2009
 */

Configuration::Configuration()
    :commandLine("Not Given")
{
  oMethod = CLI::SLO ;
  statistics_only = false;
  is_cseg = false;
  extend = false;
  verify = false;
  inputfile_given = false;
  io_at_end = true;

  inputfile = NULL;
  partitionFile = NULL;
}
/* Configuration() ENDS*/

/**
 * @Description : Desctructor
 * @Creation Date: 05.08.2009
 */
Configuration::~Configuration()
{

}
/* ~Configuration() ENDS*/


