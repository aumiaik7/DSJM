#include "NNZTag.hh"
#include <iostream>

using namespace std;


/**
 * @Description : Constructor
 * @Creation Date: 09.08.2009
 */
NNZTag::NNZTag(int _M , int _N)
  :M(_M)
  ,N(_N)
{

}
/* NNZTag() ENDS*/

/**
 * @Description : Destructor
 * @Creation Date: 09.08.2009
 */
NNZTag::~NNZTag()
{
    // cout << "Destructor for NNZTag Called" <<endl ;
}
/* ~NNZTag() ENDS*/

/**
 * @Description : Getter Function
 * @Creation Date: 09.08.2009
 */
int
NNZTag::getValue(int row, int column)
{
  return nnztag[getKey(row,column)];
}
/* getValue() ENDS*/

/**
 * @Description : Setter Method
 * @Creation Date: 09.08.2009
 */
void
NNZTag::setValue(int row, int column , int value)
{
  nnztag[getKey(row,column)] = value;
}
/* setValue() ENDS*/

/**
 * @Description : Finds the unique key combining row and column.
 * @Creation Date: 09.08.2009
 */
int
NNZTag::getKey(int row, int column)
{
  return (row -1) * N + column;
}
/* getKey() ENDS*/


