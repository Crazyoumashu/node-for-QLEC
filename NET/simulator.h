#ifndef SIMULATOR_H_INCLUDED
#define SIMULATOR_H_INCLUDED

#include <iostream>
#include <cstdlib>
#include <thread>
#include <ctime>
#include <cmath>
#include <fstream>
#include <pthread.h>
#include <unistd.h>
#include "node.h"
#include "parameters.h"

using namespace std;

class Node;
struct Center;

//std::thread thread_pool[NODE_NUMBER];

class Simulator
{
private:
    int nnode_;
    int nclusters_;
    int dropPacketNum;

    double randExp(double lambda)
    {
        double pV = 0.0;
        while(true)
        {
            pV = (double)rand()/(double)RAND_MAX;
            if (pV != 1)
            {
                break;
            }
        }
        pV = (-1.0/lambda)*log(1-pV);
        return pV;
    }

public:
    thread thread_pool[NODE_NUMBER];
    thread thread_river[NODE_NUMBER];
    static Simulator* instance_;
    Node** nodes_;
    Node* bs_;
    double* P_;
    int* N_;
    int* candidates_;
    int* clusters_heads_;
    double avg_dis_to_bs;
    double sum_energy;//init
    //bool Accept[NODE_NUMBER][NODE_NUMBER];
    int acc_num[NODE_NUMBER][NODE_NUMBER];
    int send_num[NODE_NUMBER][NODE_NUMBER];
    //FOR Q-learning
    double Q_ij[NODE_NUMBER][NODE_NUMBER];
    double P_ij[NODE_NUMBER][NODE_NUMBER];
    double R_ij[NODE_NUMBER][NODE_NUMBER];
    double V_ij[NODE_NUMBER];
    double dis_[NODE_NUMBER][NODE_NUMBER];
    int totalPacketNum;
    int acceptPacketNum;
    double totalLatency;
    bool* is_alive;
    double energy;//current consume
    //pthread_mutex_t mutex_Accept;
    pthread_mutex_t mutex_acc_num;
    pthread_mutex_t mutex_send_num;
    pthread_mutex_t mutex_write;
    pthread_mutex_t mutex_accPacket;
    pthread_mutex_t mutex_energy;
    pthread_mutex_t mutex_Q_ij;
    pthread_mutex_t mutex_R_ij;
    pthread_mutex_t mutex_V_ij;
    pthread_mutex_t mutex_latency;
    ifstream fin;
    ofstream fout;
    time_t start_time;
    time_t begin_time;
    bool thirty_percent;


    Simulator();
    ~Simulator();

    void init();
    static Simulator& instance()
    {
        return (*instance_);        // general access to scheduler
    }
    Node* node(int i)
    {
        if(i < nnode_)
            return (nodes_[i]);
        else
            return NULL;

    }
    Node** nodePtr(int i)
    {
        if(i < nnode_)
            return (nodes_+i);
        else return NULL;

    }
    Node* baseStation()
    {
        return bs_;
    }
    int nClusters()
    {
        return nclusters_;
    }
    void producePacket();
    void killNode(Node*);
    void runCHNode(int i);
    void runTCHNode(int i, int lay);
    void addEnergy();
    void k_means(int n);
    void deec(int r);
    void FCM(int n, int layer);
    void lifespan();
    void run();
};

#endif // SIMULATOR_H_INCLUDED
