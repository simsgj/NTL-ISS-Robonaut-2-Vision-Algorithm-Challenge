#include "mergepair.h"

MergePair::MergePair()
{
}


bool MergePair::operator<(const MergePair& with)
{
    return this->score < with.score;
}
