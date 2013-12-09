#ifndef MERGEPAIR_H
#define MERGEPAIR_H

#include <vector>
using namespace std;

class MergePair {

public:
  MergePair();
  bool operator<(const MergePair& with);

  size_t first;
  size_t second;
  float score;
};


#endif // MERGEPAIR_H
