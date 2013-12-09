#ifndef EDGE_H
#define EDGE_H

class Edge
{
public:
    Edge();

    bool operator< (const Edge& param) const;

    float weight;
    int from;
    int to;
};

#endif // EDGE_H
