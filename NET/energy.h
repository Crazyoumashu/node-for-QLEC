#ifndef ENERGY_H_INCLUDED
#define ENERGY_H_INCLUDED

#include "node.h"

struct Node;

//typedef enum{TX, RX} EnergyConsumption_t;
class Energy {
public:
    static void spend(Node*, Node*);
private:
    static const double energyPerBitTx;
    static const double energyPerBitRx;
};

#endif // ENERGY_H_INCLUDED
