#include "node.h"
#include "simulator.h"
#include "energy.h"


int Node::IDs_ = 0;
// Distance between two nodes
double Node :: distance(Node *n1, Node *n2) {
    double x1 = (double)(n1->location_.x);
    double y1 = (double)(n1->location_.y);
    double z1 = (double)(n1->location_.z);
    double x2 = (double)(n2->location_.x);
    double y2 = (double)(n2->location_.y);
    double z2 = (double)(n2->location_.z);
    return (sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) +pow((z2-z1), 2)));
}

int Node::reachedThreshold(){
    if(energy_ < energy_threshold_1_)
        return -1;
    else if(energy_ > energy_threshold_2_)
        return 1;
    else
        return 0;
}

int Node::enqueuePkt(Packet *p){
    if(pktQueue_.size() == queueLimit_){
        //std::cout<<"[ERROR]: Packet que is full"<<std::endl;
        return -1;
    }
    pktQueue_.push_back(*p);
    return 0;
}

int Node::dequeuePkt(){
    if(pktQueue_.size() == 0){
        //std::cout<<"[ERROR]: Packet que is empty"<<std::endl;
        return -1;
    }
    pktQueue_.erase(pktQueue_.begin());
    return 0;
}

// Send the packet p to the node n
int Node :: send(Node *n, Packet *p) {
    //cout<<p->type<<' '<<p->sourceID<<' '<<p->destID<<endl;
    Simulator &sim = Simulator::instance();

    pthread_mutex_lock(&sim.mutex_send_num);
    sim.send_num[id_][n->id_] += 1;
    pthread_mutex_unlock(&sim.mutex_send_num);
    //sim.P_ij[id_][n->id_] = (1.0 * sim.acc_num[id_][n->id_]) / sim.send_num[id_][n->id_];

    //p->forwardID = id_;
    //double trans_time = (p->size)/BIT_RATE;

    Energy::spend(this, n); //energy consumption by transmitter
    n->recv(this,p);

    return 0;
}

int Node :: recv(Node*n, Packet *p) {
    Simulator &sim = Simulator::instance();
    double ran_num = (double)(rand()/(double)RAND_MAX);
    if(this==sim.instance().bs_)
    {
        if(ran_num <=DROP_RATE)
        {
            pthread_mutex_lock(&sim.mutex_accPacket);
            sim.acceptPacketNum += 1;
            pthread_mutex_unlock(&sim.mutex_accPacket);
            pthread_mutex_lock(&sim.mutex_acc_num);
            sim.acc_num[n->id_][id_] += 1;
            pthread_mutex_unlock(&sim.mutex_acc_num);
            pthread_mutex_lock(&sim.mutex_latency);
            sim.totalLatency += time(NULL) - p->birth;
            pthread_mutex_unlock(&sim.mutex_latency);
            //sim.P_ij[n->id_][id_] = (1.0 * sim.acc_num[n->id_][id_]) / sim.send_num[n->id_][id_];
        }
        return -1;
    }
    else if(pktQueue_.size() < queueLimit_)
    {
        if(ran_num <=DROP_RATE)
        {
            Packet *p1 = new Packet;
            p1->size = COMPRESS_RATIO*PACKET_SIZE;
            p1->birth = p->birth;
            pktQueue_.push_back(*p1);
            pthread_mutex_lock(&sim.mutex_acc_num);
            sim.acc_num[n->id_][id_] += 1;
            pthread_mutex_unlock(&sim.mutex_acc_num);
            //sim.P_ij[n->id_][id_] = (1.0 * sim.acc_num[n->id_][id_]) / sim.send_num[n->id_][id_];
        }
        return 1;
    }
    else return 0;
}

void Node::forwardData_for_baseline()
{
    Simulator& sim = Simulator::instance();
    if(this->nodeType()!=TCH)
    {
        if(pktQueue_.size()!=0)
        {
            Packet *p = new Packet;
            p->size = PACKET_SIZE;
            p->birth = pktQueue_.front().birth;
            if(sim.is_alive[this->cluster()]==true)
            send(sim.nodes_[this->cluster()],p);
            else send(sim.instance().bs_, p);
            dequeuePkt();
        }
    }
    else
    {
        if(this->layers() ==1)
        {
            if(pktQueue_.size()!=0)
         {
             Packet *p = new Packet;
             p->size = PACKET_SIZE*COMPRESS_RATIO;
             p->birth = pktQueue_.front().birth;
             send(sim.instance().bs_, p);
             dequeuePkt();
         }
        }
        else
        {
            double d_thres;
            d_thres = sqrt(EPISILON_FS/EPISILON_MP);
            double min_tmp = 1000000000.0;
            double tmp;
            int ind = 0;
            for(int i = 0; i<NODE_NUMBER-1; i++)
            {
                if(sim.nodes_[i]->layers()==(this->layers()-1)&&sim.nodes_[i]->nodeType()==TCH)
                {
                    if(sim.dis_[i][id_]<d_thres)
                    tmp = EPISILON_FS*pow(sim.dis_[i][id_], 2)/sim.nodes_[i]->energy();
                    else tmp = EPISILON_MP*pow(sim.dis_[i][id_], 4)/sim.nodes_[i]->energy();

                    if(tmp < min_tmp)
                    {
                        min_tmp = tmp;
                        ind = i;
                    }
                }
            }
            if(pktQueue_.size()!=0)
         {
             Packet *p = new Packet;
             p->size = PACKET_SIZE*COMPRESS_RATIO;
             p->birth = pktQueue_.front().birth;
             send(sim.nodes_[ind], p);
             dequeuePkt();
         }
        }
    }
}

void Node::forwardData_for_kmeans()
{
    Simulator& sim= Simulator::instance();
    if(this->nodeType()!=TCH)
    {
        if(pktQueue_.size()!=0)
        {
            Packet *p = new Packet;
            p->size = PACKET_SIZE;
            p->birth = pktQueue_.front().birth;
            if(sim.is_alive[this->cluster()]==true)
            send(sim.nodes_[this->cluster()],p);
            else send(sim.instance().bs_, p);
            dequeuePkt();
        }
    }
    else
    {
         if(pktQueue_.size()!=0)
         {
             Packet *p = new Packet;
             p->size = PACKET_SIZE*COMPRESS_RATIO;
             p->birth = pktQueue_.front().birth;
             send(sim.instance().bs_, p);
             dequeuePkt();
         }
    }
}
void Node :: forwardData() {
    Simulator& sim= Simulator::instance();
    double R_t;
    if(this->nodeType()!=TCH)
    {
        if(pktQueue_.size()!=0)
        {
            double Max = -100000;
            int ch_id;
            Packet *p = new Packet;//(pktQueue_.front());
            p->size = PACKET_SIZE;
            p->birth = pktQueue_.front().birth;
            for(int i = 0; i<NODE_NUMBER; i++)
            {
                if(sim.nodes_[i]->nodeType()== TCH || sim.nodes_[i]->nodeType() == BS)
                {
                    if(sim.nodes_[i]->nodeType()== TCH)
                    {
                        pthread_mutex_lock(&sim.mutex_R_ij);
                        sim.R_ij[id_][i] = -G + ALPHA1*(this->energy()+sim.nodes_[i]->energy())-ALPHA2*sim.dis_[id_][i];
                        pthread_mutex_unlock(&sim.mutex_R_ij);
                    }
                    else
                    {
                        pthread_mutex_lock(&sim.mutex_R_ij);
                        sim.R_ij[id_][i] = -G + ALPHA1*(this->energy()+sim.nodes_[i]->energy())-ALPHA2*sim.dis_[id_][i] - LL;
                        pthread_mutex_unlock(&sim.mutex_R_ij);
                    }

                    pthread_mutex_lock(&sim.mutex_R_ij);
                    sim.R_ij[id_][id_] = -G + Beta1*this->energy()-Beta2*sim.dis_[id_][i];
                    pthread_mutex_unlock(&sim.mutex_R_ij);

                    R_t = sim.P_ij[id_][i]*sim.R_ij[id_][i]+(1-sim.P_ij[id_][i])*sim.R_ij[id_][id_];

                    pthread_mutex_lock(&sim.mutex_Q_ij);
                    sim.Q_ij[id_][i] = R_t + DISCOUNT_FACTOR*(sim.P_ij[id_][i]*sim.V_ij[i]+(1-sim.P_ij[id_][i])*sim.V_ij[id_]);
                    pthread_mutex_unlock(&sim.mutex_Q_ij);
                    if(sim.Q_ij[id_][i]>Max)
                    {
                        Max = sim.Q_ij[id_][i];
                        ch_id = i;
                    }
                }
            }
            pthread_mutex_lock(&sim.mutex_V_ij);
            sim.V_ij[id_] = Max;
            pthread_mutex_unlock(&sim.mutex_V_ij);
            this->cluster(ch_id);
            send(sim.nodes_[this->cluster()],p);
            //pktQueue_.erase(pktQueue_.begin());
            dequeuePkt();

        }

    }
    else
    {
        if(pktQueue_.size()!=0)
        {
            Packet *p = new Packet;//(pktQueue_.front());
            p->size = PACKET_SIZE*COMPRESS_RATIO;
            p->birth = pktQueue_.front().birth;
            pthread_mutex_lock(&sim.mutex_R_ij);
            sim.R_ij[id_][NODE_NUMBER-1] = -G + ALPHA1*(this->energy()+sim.nodes_[NODE_NUMBER-1]->energy())-ALPHA2*sim.dis_[id_][NODE_NUMBER-1];
            pthread_mutex_unlock(&sim.mutex_R_ij);

            pthread_mutex_lock(&sim.mutex_R_ij);
            sim.R_ij[id_][id_] = -G + Beta1*this->energy()-Beta2*sim.dis_[id_][NODE_NUMBER];
            pthread_mutex_unlock(&sim.mutex_R_ij);

            R_t = sim.P_ij[id_][NODE_NUMBER]*sim.R_ij[id_][NODE_NUMBER]+(1-sim.P_ij[id_][NODE_NUMBER])*sim.R_ij[id_][id_];

            pthread_mutex_lock(&sim.mutex_Q_ij);
            sim.Q_ij[id_][NODE_NUMBER-1] = R_t + DISCOUNT_FACTOR*(sim.P_ij[id_][NODE_NUMBER-1]*sim.V_ij[NODE_NUMBER-1]+(1-sim.P_ij[id_][NODE_NUMBER-1])*sim.V_ij[id_]);
            pthread_mutex_unlock(&sim.mutex_Q_ij);

            pthread_mutex_lock(&sim.mutex_V_ij);
            sim.V_ij[id_] = sim.Q_ij[id_][NODE_NUMBER-1];
            pthread_mutex_unlock(&sim.mutex_V_ij);

            send(sim.instance().bs_, p);
            //pktQueue_.erase(pktQueue_.begin());
            dequeuePkt();
        }
    }

}
