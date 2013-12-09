#include "edge.h"

Edge::Edge()
{
}

bool Edge::operator <(const Edge& param) const {
    return this->weight < param.weight;
}
