#ifndef NODE_H_INCLUDED
#define NODE_H_INCLUDED
#include <vector>
#include <cstring>
#include "net.h"
#include "mac.h"
#include <math.h>
#include <stdio.h>


// Node type
typedef enum{CH, BS, TCH} NodeType_t;

struct Center
{
    double x;
    double y;
    double z;
};
// Represents a wireless sensor node
class Node{
public:
    std::vector<Packet> pktQueue_;

    Node(){
        id_ = IDs_++;
        nodeType_ = CH;
        maxEnergy_ = INIT_ENERGY;
        energy_ = INIT_ENERGY;
        location_.x = location_.y = location_.z = 0;
        not_cluster_times_ = MAX_T;
        queueLimit_ = DEFAULT_QUEUE_LIMIT;
        //memset(Vfuction, 0, sizeof(Vfuction));
        //transmissionRange_ = MAX_RANGE;
        pktQueue_.clear();
    }
    //void init(int nneighbor);

    // Accessors
    int id() { return id_;}
    double maxEnergy() { return maxEnergy_;}
    void maxEnergy(double val) {maxEnergy_ = val;}
    double lastEnergy() {return lastEnergy_;}
    void lastEnergy(double val) {lastEnergy_ = val;}
    double energy() { return energy_;}
    void energy(double val) {energy_ = val;}
    //double transmissionRange() { return transmissionRange_; }
    //void transmissionRange(double val) { transmissionRange_ = val; }
    Cordinate location() { return location_; }
    void location(double x, double y, double z) { location_.x = x; location_.y = y; location_.z = z; }
    //Node* nextHop() { return nextHop_; }
    void nodeType(NodeType_t t) { nodeType_ = t; }
    NodeType_t nodeType() { return nodeType_; }
    int cluster() { return clusterId_; }
    void cluster(int cl) { clusterId_ = cl; }
    int not_cluster_times() {return not_cluster_times_;}
    void not_cluster_times(int times) {not_cluster_times_ = times;}
    int layers(){return layers_;}
    void layers(int layer_num){layers_ = layer_num;}


    //Functions
    //void nextHop(Node *n) {nextHop_ = n;}
    //void eventData(const char *data);
    double distance(Node*, Node*); // Distance between two nodes
    //Node* clusterHead() { return neighbours_[clusterId_]; }
    int reachedThreshold(); // Check if energy level has reached a threshold
    //void selectNextHop(int k);
    void forwardData(); // This calls selectNextHop() and forwards packet to it.
    void forwardData_for_kmeans();
    void forwardData_for_baseline();

    int enqueuePkt(Packet*);
    int dequeuePkt();
    int send(Node*, Packet*);    // Send the packet p to the node n(for special purposes)
    int recv(Node*, Packet*);    // Recieve and enqueu packet p
    //int broadcast(Packet*);    // Broadcast packet p
    //int notifyRelax();    // Broadcast a relaxation packet

private:
    static int IDs_;  // ID generated and ID count
    double maxEnergy_;
    int sendPacket; // Total packet send from i
    //double Vfuction[K+2]; // Q-learning

    //MAC layer properties
    Cordinate location_;
    double energy_;
    double lastEnergy_;
    //double harvest_;
    double energy_threshold_1_ = ENERGY_THRESHOLD_1;
    double energy_threshold_2_ = ENERGY_THRESHOLD_2;
    //double transmissionRange_;
    unsigned int queueLimit_; //Max pktqueue size after which the packets will be dropped.

    //Net layer properties
    int id_; //ID or address
    NodeType_t nodeType_; //CH,ACH,NCH,BS,TCH
    int clusterId_; // Cluster number (NOTE: ClusterID == ID of CH node)
    int not_cluster_times_;
    int layers_;
};

#endif // NODE_H_INCLUDED
