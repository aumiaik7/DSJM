#ifndef NNZTAG_HH
#define NNZTAG_HH
#include <map>

using namespace std;

struct ltint
{
  bool operator() (const int i1, const int i2) const
  {
    return (i1 < i2);
  }
};
class NNZTag
{
private:
  int getKey(int row, int column);
  int M, N;
  map<int, int, ltint> nnztag;
public:
  NNZTag(int m, int n);
  virtual ~NNZTag();

  int getValue(int row,int column);
  void setValue(int row,int column, int value);
};

#endif
