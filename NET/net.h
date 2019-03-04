#ifndef NET_H_INCLUDED
#define NET_H_INCLUDED
#include <cstring>
#include "node.h"
#include "parameters.h"
#include <iostream>


struct Packet {
    int size;
    time_t birth;
};

#endif // NET_H_INCLUDED
