#include "simulator.h"

int main(int argc, const char * argv[])
{
    (new Simulator);
    Simulator &sim = Simulator::instance();
    sim.init();
    sim.run();

    /*double total_energy[51][134];
    double consume_energy[52][134];
    double consume_ratio[52][134];
    memset(total_energy, 0, sizeof(total_energy));
    memset(consume_energy, 0, sizeof(consume_energy));
    memset(consume_ratio, 0, sizeof(consume_ratio));*/
    //cout<<100<<endl;

    /*for(int i = 0; i<NODE_NUMBER; i++)
    {
        for(int j = 0; j<NODE_NUMBER; j++)
        {
            if(j == NODE_NUMBER-1)
                cout<<sim.P_ij[i][j]<<endl;
            else cout<<sim.P_ij[i][j]<< " ";
        }

    }*/
    if (sim.totalPacketNum > 0)
    {
        cout << "Packet delivery ratio: " << sim.acceptPacketNum / double(sim.totalPacketNum) << endl;
        sim.fout<<sim.acceptPacketNum / double(sim.totalPacketNum)<<endl;
        //cout << "Latency per Packet: "<< sim.totalLatency /sim.acceptPacketNum<<endl;
        //cout<< "Energy Cost: "<<sim.consume_energy<<endl;
    }
    else
    {
        cout << "Total packet number is 0!!!" << endl;
    }
    /*int tmp_x,tmp_y;
    for(int i = 0; i< NODE_NUMBER-1; i++)
    {
        tmp_x = (int)sim.nodes_[i]->location().x;
        tmp_y = (int)sim.nodes_[i]->location().y;
        //cout<<tmp_x<<" "<<tmp_y<<" "<<i<<endl;
        total_energy[tmp_x][tmp_y] += sim.nodes_[i]->maxEnergy();
        consume_energy[tmp_x][tmp_y] += sim.nodes_[i]->maxEnergy() - sim.nodes_[i]->energy();
    }


    for(int i = 18; i<51; i++)
    {
        for(int j = 75; j<133; j++)
        {
            if(total_energy[i][j]!=0)
            {
                consume_ratio[i][j] = consume_energy[i][j]/total_energy[i][j];
            }
            sim.fout<<j<<","<<i<<","<<consume_ratio[i][j]<<endl;
        }
    }*/


    /*for (int i = 0; i < NODE_NUMBER-1; ++i)
    {
        cout << "Rest energy of node " << i << ": " << sim.nodes_[i]->energy() << endl;
    }*/

    return 0;
}
