#ifndef QUADRUPLEINDEXES_H
#define QUADRUPLEINDEXES_H

class QuadrupleIndexes
{
public:
  QuadrupleIndexes();

  bool operator<(const QuadrupleIndexes& with) const;

  unsigned long i;
  unsigned long j;
  unsigned long k;
  unsigned long l;
  float score;
};

#endif // QUADRUPLEINDEXES_H
