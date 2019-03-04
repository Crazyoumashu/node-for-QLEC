#include "simulator.h"
#define min(a,b) (a>b)?b:a

Simulator* Simulator :: instance_;

Simulator::Simulator()
{
    nnode_ = NODE_NUMBER;
    avg_dis_to_bs = 0.0;
    sum_energy = 0;
    energy = 0;
    totalPacketNum = 0;
    acceptPacketNum = 0;
    totalLatency = 0.0;
    thirty_percent = false;
    instance_ = this;
    nclusters_ = 0;
    P_ = new double[nnode_-1];
    N_ = new int[nnode_-1];
    candidates_ = new int[nnode_-1];
    clusters_heads_ = new int[nnode_-1];
    nodes_ = new Node*[nnode_];
    is_alive = new bool[nnode_];

    //memset(Accept, 0, sizeof(Accept));
    memset(acc_num, 0, sizeof(acc_num));
    memset(send_num, 0, sizeof(send_num));
    memset(Q_ij, 0, sizeof(Q_ij));
    memset(R_ij, 0, sizeof(R_ij));
    memset(P_ij, 0, sizeof(P_ij));
    memset(V_ij, 0, sizeof(V_ij));

    //pthread_mutex_init(&mutex_Accept, NULL);
    pthread_mutex_init(&mutex_acc_num, NULL);
    pthread_mutex_init(&mutex_send_num, NULL);
    pthread_mutex_init(&mutex_write, NULL);
    pthread_mutex_init(&mutex_accPacket, NULL);
    pthread_mutex_init(&mutex_energy, NULL);
    pthread_mutex_init(&mutex_Q_ij, NULL);
    pthread_mutex_init(&mutex_R_ij, NULL);
    pthread_mutex_init(&mutex_V_ij, NULL);
    pthread_mutex_init(&mutex_latency, NULL);

    fout.open("out.txt");
    //fin.open("loc.txt",ios::in);
}

Simulator::~Simulator()
{
    for (int i = 0; i < nnode_; ++i)
    {
        delete nodes_[i];
        nodes_[i] = NULL;
    }
    delete []nodes_;
    delete []is_alive;
}

void Simulator :: init()
{
    double x, y, z;
    //double d_x, d_y, d_z, d_energy;
    int *init_energy;
    //int sum = 0;
    init_energy = new int[nnode_];
    for (int i = 0; i < nnode_-1; i++)
    {
        init_energy[i] = 5;
        sum_energy += init_energy[i];
    }
    //init_energy[nnode_] = 5000*nnode_ - sum;
    for (int i = 0; i < nnode_; ++i)
    {
        nodes_[i] = new Node;
        is_alive[i] = 1;
        if (i < nnode_-1)
        {
            nodes_[i]->nodeType(CH);

            nodes_[i]->not_cluster_times(MIN_DISTANCE);
            double current_r = rand()*1.0/RAND_MAX*R;
            int current_ag1 = rand()%360;
            int current_ag2 = rand()%360;
            //x = rand() % (MAX_X+1);
            //y = rand() % (MAX_Y+1);
            //z = rand() % (MAX_Z-MIN_Z+1)+MIN_Z;
            z = sin(current_ag1)*current_r;
            x = cos(current_ag1)*cos(current_ag2)*current_r;
            y = cos(current_ag1)*sin(current_ag2)*current_r;
            //fin>>d_x>>d_y>>d_z>>d_energy;
            //cout<<d_x<<" "<<d_y<<" "<<d_z<<endl;
            /*x = (double)d_x;
            y = (double)d_y;
            z = (double)d_z;
            if(i == 2886)
            {
                int tmp = x;
                x = y;
                y = tmp;
            }*/
            nodes_[i]->location(x,y,z);
            nodes_[i]->energy(init_energy[i]);
            nodes_[i]->maxEnergy(init_energy[i]);
            nodes_[i]->lastEnergy(init_energy[i]);
            sum_energy += init_energy[i];
        }
        else
        {
            nodes_[i]->nodeType(BS);
            nodes_[i]->not_cluster_times(0);
            nodes_[i]->energy(5);
            nodes_[i]->location(0.5*(MAX_X+MIN_X), 0.5*(MAX_Y+MIN_Y), 0.5*MAX_Z);
            bs_ = nodes_[i];
        }
    }
    for(int i = 0;i <nnode_; i++)
    {
        for(int j = 0; j<nnode_; j++)
        {
            double x1 = (double)(nodes_[i]->location().x);
            double y1 = (double)(nodes_[i]->location().y);
            double z1 = (double)(nodes_[i]->location().z);
            double x2 = (double)(nodes_[j]->location().x);
            double y2 = (double)(nodes_[j]->location().y);
            double z2 = (double)(nodes_[j]->location().z);
            dis_[i][j] = sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) +pow((z2-z1), 2));
        }
    }
    for(int i = 0; i < nnode_-1; ++i)
    {
        avg_dis_to_bs += dis_[i][nnode_-1];
        //cout<<dis_[i][nnode_-1]<<endl;
        //cout<<avg_dis_to_bs<<endl;
    }
    avg_dis_to_bs /= nnode_-1;

    for(int i = 0; i<nnode_-1; i++)
    {
        double x1 = (double)(nodes_[i]->location().x);
        double y1 = (double)(nodes_[i]->location().y);
        double z1 = (double)(nodes_[i]->location().z);
        double tmp = sqrt(pow(x1, 2) + pow(y1, 2) + pow(z1, 2));
        int lay = (int)(tmp/(R/LAYER_NUM));
        if(lay!=LAYER_NUM)
        {
            nodes_[i]->layers(lay+1);
        }
        else nodes_[i]->layers(lay);
    }
    //cout<<avg_dis_to_bs<<endl;

    for(int i = 0; i<NODE_NUMBER; i++)
    {
        for(int j = 0; j<NODE_NUMBER; j++)
        {
            P_ij[i][j] = 0.99;
        }
    }
}

void Simulator :: FCM(int n, int layer)
{
    std::vector<Center>center(n);
    Node** node_g;
    double** U_ij;
    double** dis;
    bool** be_cluster;
    be_cluster = new bool*[n];
    dis = new double*[n];
    double sum_x, sum_y, sum_z;
    double sum_u;
    double val;
    double val_last;
    val_last = -20;
    int cnt = 0;
    for(int i = 0; i<nnode_ -1; i++)
    {
        if(nodes_[i] ->layers()==layer)
        {
            cnt++;
        }
    }
    node_g = new Node*[cnt];
    U_ij = new double*[n];
    for(int i = 0; i<n; i++)
    {
        U_ij[i] = new double[cnt];
        dis[i] = new double[cnt];
        be_cluster[i] = new bool[cnt];
    }
    int index = 0;
    for(int i = 0; i<nnode_ - 1; i++)
    {
        if(nodes_[i] ->layers()==layer)
        {
            node_g[index] = nodes_[i];
            index++;
        }
    }
    for(int i = 0; i<n; i++)
    {
        for(int j = 0; j<cnt; j++)
        {
            if(i == 0)
                U_ij[i][j] = 0.4;
            else U_ij[i][j] = 0.6;

            be_cluster[i][j] = 0;
        }
    }
    double delta = 10000;
    for(int i = 0; i<n; i++)
    {
        center[i].x = 0;
        center[i].y = 0;
        center[i].z = 0;
    }
    while(delta > 10)
    {
        for(int i = 0; i<n; i++)
        {
            sum_x = 0;
            sum_y = 0;
            sum_z = 0;
            for(int j = 0; j<cnt; j++)
            {
                sum_x += U_ij[i][j]*node_g[j]->location().x;
                sum_y += U_ij[i][j]*node_g[j]->location().y;
                sum_z += U_ij[i][j]*node_g[j]->location().z;
            }
            center[i].x = sum_x/cnt;
            center[i].y = sum_y/cnt;
            center[i].z = sum_z/cnt;
        }
        for(int i = 0; i<n; i++)
        {
            double x1 = center[i].x;
            double y1 = center[i].y;
            double z1 = center[i].z;
            for(int j = 0; j<cnt; j++)
            {
                double x2 = node_g[j]->location().x;
                double y2 = node_g[j]->location().y;
                double z2 = node_g[j]->location().z;
                dis[i][j] = sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) +pow((z2-z1), 2));
            }
        }

        for(int i = 0; i<n; i++)
        {
            for(int j = 0; j< cnt; j++)
            {
                sum_u = 0;
                for(int k =0; k<n; k++)
                {
                    sum_u += pow(dis[i][j]/dis[k][j], 2.0/(M-1));
                }
                U_ij[i][j] = 1.0/sum_u;
            }
        }
        val = 0;
        for(int i = 0; i<n; i++)
        {
            for(int j = 0; j<cnt; j++)
            {
                val += pow(U_ij[i][j], M)*pow(dis[i][j], 2);
            }
        }
        delta = val - val_last;
        val_last = val;
    }
    double MAX;
    int ind;
    for(int j = 0; j<cnt; j++)
    {
        MAX = 0;
        ind = 0;
        for(int i = 0; i<n; i++)
        {
            if(U_ij[i][j]>MAX)
            {
                MAX = U_ij[i][j];
                ind = i;
            }
        }
        be_cluster[ind][j] = 1;
    }
    for(int i = 0; i<n; i++)
    {
        double MIN = MIN_DISTANCE;
        ind = 0;
        for(int j =0; j<cnt; j++)
        {
            if(be_cluster[i][j]==1)
            {
                if(dis[i][j]< MIN)
                {
                    ind = j;
                    MIN = dis[i][j];
                }
            }
        }
        nodes_[node_g[ind]->id()]->nodeType(TCH);
        for(int k = 0; k<cnt; k++)
        {
            if(be_cluster[i][k]==1&&k!=ind)
            {
                nodes_[node_g[k]->id()]->cluster(node_g[ind]->id());
                nodes_[node_g[k]->id()]->nodeType(CH);
            }
        }
    }

}

void Simulator :: k_means(int n)
{
    std::vector<Center>center(n);
    std::vector<Center>last_center(n);
    double delta = 1000;
    double**dis;
    bool**becluster;
    bool**becluster_last;
    int *cnt;
    int *min_dis;
    int sum_x, sum_y, sum_z;
    cnt = new int[n];
    min_dis = new int[nnode_-1];
    dis = new double*[nnode_-1];
    becluster = new bool*[nnode_-1];
    becluster_last = new bool*[nnode_-1];
    for(int i = 0; i<nnode_-1; i++)
    {
        dis[i] = new double[n];
        becluster[i] = new bool[n];
        becluster_last[i] = new bool[n];
    }
    for(int i = 0; i<nnode_-1; i++)
    {
        for(int j = 0; j<n; j++)
        {
            becluster[i][j] = 0;
            becluster_last[i][j] = 0;
        }
    }
    //init centers
    for(int i = 0; i< n; i++)
    {
        center[i].x = rand() % (MAX_X+1);
        center[i].y = rand() % (MAX_Y+1);
        center[i].z = rand() % (MAX_Z-MIN_Z+1)+MIN_Z;
        cnt[i] = 0;
    }
    while(delta > 10)//convegence
    {
        delta = 0;
        //divide node
        for(int i = 0; i<nnode_-1; i++)
        {
            for(int j = 0; j<n; j++)
            {
                double x1 = (double)(nodes_[i]->location().x);
                double y1 = (double)(nodes_[i]->location().y);
                double z1 = (double)(nodes_[i]->location().z);
                double x2 = (double)(center[j].x);
                double y2 = (double)(center[j].y);
                double z2 = (double)(center[j].z);
                dis[i][j] = sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) +pow((z2-z1), 2));
            }
        }
        int index;
        for(int i = 0; i< nnode_-1; i++)
        {
            if(is_alive[i])
            {
                min_dis[i] = MIN_DISTANCE;
                index = 0;
                for(int j = 0; j<n; j++)
                {
                    if(dis[i][j]<min_dis[i])
                    {
                        min_dis[i] = dis[i][j];
                        index = j;
                    }
                }
                    becluster[i][index] = 1;
            }
        }
        for(int i = 0; i< nnode_-1; i ++)
        {
            for(int j = 0; j<n; j++)
            {
                becluster_last[i][j] = becluster[i][j];
            }
        }
        //find new centers
        for(int j =0; j<n; j++)
        {
            sum_x = 0;
            sum_y = 0;
            sum_z = 0;
            for(int i = 0; i<nnode_-1; i++)
            {
                if(is_alive[i])
                {
                    if(becluster[i][j] == 1)
                    {
                    cnt[j]+= 1;
                    sum_x += nodes_[i]->location().x;
                    sum_y += nodes_[i]->location().y;
                    sum_z += nodes_[i]->location().z;
                    becluster[i][j] = 0;
                    }
                }
            }
            last_center[j].x = center[j].x;
            last_center[j].y = center[j].y;
            last_center[j].z = center[j].z;
            if(cnt[j]!=0)
            {
                center[j].x = sum_x/cnt[j];
                center[j].y = sum_y/cnt[j];
                center[j].z = sum_z/cnt[j];
            }
            cnt[j] = 0;
        }
        //calculate delta
        for(int j = 0; j<n; j++)
        {
             double x1 = (double)(last_center[j].x);
             double y1 = (double)(last_center[j].y);
             double z1 = (double)(last_center[j].z);
             double x2 = (double)(center[j].x);
             double y2 = (double)(center[j].y);
             double z2 = (double)(center[j].z);
             double tmp = sqrt( pow((x2-x1), 2) + pow((y2-y1), 2) +pow((z2-z1), 2));
             delta += tmp;
        }
    }//end
    //we choose the node that is closest to center to be the clusterhead
    for(int j = 0; j<n; j++)
    {
        int min = MIN_DISTANCE;
        int ind = 0;
        for(int i = 0; i<nnode_-1; i++)
        {
            if(is_alive[i])
            {
                if(becluster_last[i][j]==1)
                {
                if(dis[i][j]<min)
                {
                    min = dis[i][j];
                    ind = i;
                }
                }
            }
        }
        nodes_[ind]->nodeType(TCH);
        for(int k = 0; k<nnode_-1; k++)
        {
            if(is_alive[k])
            {
                if(becluster_last[k][j]==1)
                {
                if(k!=ind)
                {
                    nodes_[k]->cluster(nodes_[ind]->id());
                    nodes_[k]->nodeType(CH);
                }

                }
            }
        }
    }
}
void Simulator :: deec(int r)
{
    int final_count = 0;
    double k_opt, p_opt, dis_c;
    int N = nnode_-1;
    //to be modified
    //k_opt = (sqrt(N)/sqrt(2*PI))*sqrt(EPISILON_FS/EPISILON_MP)*(M1/pow(avg_dis_to_bs,2));
    double tmpt1;//tmpt2;
    //tmpt2 = (r*1.0)/ROUNDS;
    tmpt1 = (8*PI*N*EPISILON_FS)/(15*EPISILON_MP);
    k_opt = 3.0/(4*PI)*pow(tmpt1, 3.0/5)*pow(M1, 6.0/5)/pow(avg_dis_to_bs, 12.0/5);
    p_opt = k_opt/N;
    dis_c = M1*pow(3/(4*PI*k_opt), 1.0/3);
    //cout<<avg_dis_to_bs<<endl;
    cout<<k_opt<<" "<<p_opt<<" "<<dis_c<<endl;
    double avg_energy = sum_energy/N * (1-((r*1.0)-1.0)/ROUNDS);
    for(int i = 0; i<N; i++)
    {
        //cout<<nodes_[i]->energy()<<" "<<avg_energy<<endl;
        P_[i] = p_opt*(nodes_[i]->energy()/avg_energy);
        //cout<<P_[i]<<endl;
        N_[i] = (int)(1/P_[i]);
        if(N_[i]<=1)
            N_[i] = 2;
        //cout<<N_[i]<<endl;
    }
    //cout<<tmpt2<<endl;
    int count = 0;
    int tmp;
    //bool *flag1;//candidate
    //bool *flag2;//clusterhead
    bool *flag3;//candidate but not a cluster
    flag3 = new bool[N];
    //flag2 = new int[N];
    //flag1 = new int[N];
    for(int i = 0; i<N; i++)
    {
        //flag1[i] = 0;
        //flag2[i] = 0;
        flag3[i] = 0;
    }
    for(int i = 0; i<N; i++)
    {
        if(nodes_[i]->not_cluster_times()>= N_[i]&&nodes_[i]->energy()>=(1-pow(r*1.0/ROUNDS,2))*nodes_[i]->maxEnergy()&&is_alive[i]==1)
        {
            candidates_[count] = i;
            count++;
        }
        else
        {
            tmp = nodes_[i]->not_cluster_times()+1;
            nodes_[i]->not_cluster_times(tmp);
            nodes_[i]->nodeType(CH);
        }
    }
    //cout<<count<<endl;
    int cluster_count = 0;
    double Threshold;
    for(int i = 0; i<count; i++)
    {
        double ran = ((double)(rand()))/RAND_MAX;
        Threshold = P_[candidates_[i]]/(1-P_[candidates_[i]]*(r%N_[candidates_[i]]));
        //cout<<Threshold<<endl;
        if(ran < Threshold)
        {
            clusters_heads_[cluster_count] = candidates_[i];
            cluster_count += 1;
        }
        else
        {
            tmp = nodes_[candidates_[i]]->not_cluster_times()+1;
            nodes_[candidates_[i]]->not_cluster_times(tmp);
        }
    }
    for(int i = 0; i< cluster_count; i++)
    {
        if(i!= cluster_count -1)
        {
            if(flag3[clusters_heads_[i]]==0)
            {
                for(int j = i+1; j<count; j++)
            {

                if(dis_[clusters_heads_[i]][clusters_heads_[j]]<=dis_c)
                {
                    if(nodes_[clusters_heads_[i]]->energy()>nodes_[clusters_heads_[j]]->energy())
                    {
                        flag3[clusters_heads_[j]] = 1;
                        break;
                    }
                    //else flag3[clusters_heads_[i]] = 1;
                }
            }
            }

        }

    }
    for(int i = 0; i<cluster_count; i++)
    {
        if(flag3[clusters_heads_[i]]==1)
        {
            tmp = nodes_[clusters_heads_[i]]->not_cluster_times()+1;
            nodes_[clusters_heads_[i]]->not_cluster_times(tmp);
            nodes_[clusters_heads_[i]]->nodeType(CH);
        }
        else
        {
            nodes_[clusters_heads_[i]]->nodeType(TCH);
            nodes_[clusters_heads_[i]]->not_cluster_times(0);
            final_count++;

        }
    }
    nclusters_ = final_count;

}

void Simulator :: lifespan()
{
    time_t cur_time;
    int old_flag = -1;
    while (true)
    {
        cur_time = time(NULL);
        if(cur_time - start_time > SIMULATE_TIME*100)
        {
            return;
        }

        int cnt_dead = 0;
        int cnt_live;
    for(int i = 0; i<nnode_-1; i++)
    {
        if(is_alive[i]==0)
        {
            cnt_dead++;
        }
    }
    cnt_live = nnode_ - cnt_dead;
    if((cur_time - begin_time)%60 == 0 &&(cur_time -begin_time)/60 !=old_flag)
    {
        cout<<cnt_dead<<" dead "<< cnt_live<< " alive "<<endl;
        old_flag = (cur_time -begin_time)/60;
    }

    double dead_rate = cnt_dead/double(nnode_);
    if(dead_rate >= 0.01&&thirty_percent == false)
    {
        int life = cur_time - begin_time;
        cout<<" Lifespan: "<< life <<" Seconds "<<endl;
        fout<<" Lifespan: "<< life <<" Seconds "<<endl;
        thirty_percent = true;
    }
    }

}

void Simulator :: producePacket()
{
    //cout<<1000<<endl;
    time_t cur_time;
    int interval = 0;
    while (true)
    {
        cur_time = time(NULL);
        if (cur_time - start_time > SIMULATE_TIME * 60)
        {
            return;
        }
            int source_id = rand()%(nnode_-1);
            while(is_alive[source_id]==0)
            source_id = rand()%(nnode_-1);
            Packet *p = new Packet;
            p->size = PACKET_SIZE;
            p->birth = time(NULL);
            if(nodes_[source_id]->nodeType()==TCH)
            {
                p->size = COMPRESS_RATIO*p->size;
            }
            nodes_[source_id]->enqueuePkt(p);
            totalPacketNum++;

        //pthread_mutex_lock(&mutex_write);
        //fout << totalPacketNum <<" produced! It's source ID :"<<source_id<<endl;
        //pthread_mutex_unlock(&mutex_write);
        interval = int(randExp(1.0/LAMBDA));
        while(interval == 0)
            interval = int(randExp(1.0/LAMBDA));
        while (time(NULL)-cur_time < interval && cur_time - start_time <= SIMULATE_TIME * 60);
    }
}

void Simulator :: runCHNode(int i)
{
    //cout<<100<<endl;
    time_t cur_time;
    while (true)
    {
        cur_time = time(NULL);
        if (cur_time - start_time > SIMULATE_TIME * 60)
        {
            //cout << "CH nodes Simulate end!!!" << endl;
            return;
        }
        if (is_alive[i])
        {
            //nodes_[i]->forwardData();
            //nodes_[i]->forwardData_for_kmeans();
            nodes_[i]->forwardData_for_baseline();
        }
        else
        {
            if(nodes_[i]->reachedThreshold()==1)
                is_alive[i] = true;
        }
    }
}
void Simulator:: runTCHNode(int i, int lay)
{
    //cout<<100<<endl;

    time_t cur_time;
    while (true)
    {
        cur_time = time(NULL);
        if (cur_time - start_time > SIMULATE_TIME * (210-30*lay))
        {
            //cout << "TCH nodes Simulate end!!!" << endl;
            return;
        }
        if (is_alive[i])
        {
            //nodes_[i]->forwardData();
            //nodes_[i]->forwardData_for_kmeans();
            nodes_[i]->forwardData_for_baseline();
        }
        else
        {
            if(nodes_[i]->reachedThreshold()==1)
                is_alive[i] = true;
        }
    }
}

void Simulator :: addEnergy()
{
    //cout<<99<<endl;
    time_t cur_time;
    time_t old_time;
    old_time = start_time;
    int old_flag = -1;
    while (true)
    {
        cur_time = time(NULL);
        if (cur_time - start_time > SIMULATE_TIME * 180)
            return;
        if (cur_time - old_time >= 1)
        {
            for (int i = 0; i < nnode_-1; ++i)
            {
                if (cur_time-start_time > SIMULATE_TIME*30)
                continue;
                nodes_[i]->energy(min(INIT_ENERGY,nodes_[i]->energy()));
            }
        }

        if ((cur_time - start_time) % 20 == 0 && int((cur_time - start_time) / 20) != old_flag)
        {
            old_flag = int((cur_time - start_time) / 20);
            pthread_mutex_lock(&mutex_write);
            fout << "Time: " << (cur_time-start_time)/60.0 << "min " << "Use: " << energy << endl;
            pthread_mutex_unlock(&mutex_write);
        }
        old_time = cur_time;
    }
}

void Simulator :: run()
{
    int cnt_ch;
    int cnt_tch;
    begin_time = time(NULL);
    for(int k = 1; k<=LAYER_NUM; k++)
    {
        FCM(2, k);
    }
    for(int r = 1; r<=ROUNDS; r++)
    {
        cnt_ch = 0;
        cnt_tch = 0;
        /*for(int i = 0; i<nnode_-1; i++)
        {
            nodes_[i]->pktQueue_.clear();
        }*/
        start_time = time(NULL);
        cout <<" Round : "<<r<< " Start time: " << start_time << endl;
        if(r>1)
        {
            for(int i = 0; i<nnode_-1; i++)
            {
                if(nodes_[i]->nodeType()==TCH)
                {
                    if(nodes_[i]->energy()<= FACTOR*nodes_[i]->lastEnergy())
                    {
                        double max_tmp = 0;
                        double tmp;
                        int ind = 0;
                        for(int j = 0; j<nnode_-1; i++)
                        {
                            if(nodes_[j]->cluster()==i)
                            {
                                tmp = nodes_[j]->energy()/dis_[i][j];
                                if(tmp>max_tmp)
                                {
                                    max_tmp = tmp;
                                    ind = j;
                                }
                            }
                        }
                        nodes_[i]->nodeType(CH);
                        nodes_[i]->cluster(ind);
                        nodes_[ind]->nodeType(TCH);
                        for(int k = 0; k<nnode_-1; k++)
                        {
                            if(nodes_[k]->cluster()==i&&k!=ind)
                            {
                                nodes_[k]->cluster(ind);
                            }
                        }
                    }
                }
            }
        }
        //deec(r);
        //k_means(nclusters_);
        //cout <<"right_deec "<<endl;
        thread thread_energy(bind(&Simulator::addEnergy,this));
        thread thread_packet(bind(&Simulator::producePacket,this));
        thread thread_life(bind(&Simulator::lifespan,this));
        for(int i = 0; i<nnode_-1; i++)
        {
            if(nodes_[i]->nodeType()!= TCH)
            {
                thread_pool[cnt_ch] = thread(bind(&Simulator::runCHNode,this, i));
                cnt_ch++;
            }
        }
        //cout<<cnt_ch<<endl;

        for (int i = 0; i < cnt_ch; ++i)
        {
            thread_pool[i].join();
        }

        thread_packet.join();

        for(int j =LAYER_NUM; j>=1; j--)
        {
           for(int i = 0; i<nnode_-1; i++)
        {
            //cout<<nodes_[i]->id()<<endl;
            if(nodes_[i]->nodeType()==TCH&&nodes_[i]->layers()==j)
            {
                thread_river[cnt_tch] = thread(bind(&Simulator::runTCHNode,this, i, j));
                cnt_tch++;
            }

        }
        }

        //cout<<cnt_tch<<endl;
        for(int i = 0; i< cnt_tch; i++)
        {

            thread_river[i].join();
        }
        thread_energy.join();
        thread_life.join();
    }

    //fout.close();
    //fin.close();
}

void Simulator :: killNode(Node* n)
{
    is_alive[n->id()] = 0;
    n->pktQueue_.clear();
}
