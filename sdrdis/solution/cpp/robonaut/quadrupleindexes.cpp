#include "quadrupleindexes.h"

QuadrupleIndexes::QuadrupleIndexes()
{
}

bool QuadrupleIndexes::operator<(const QuadrupleIndexes& with) const
{
    return this->score < with.score;
}
