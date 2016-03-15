#include <Two_Strain_Percolation_Sim.h>

int main() {
    
    // Construct Network
    Network net("name", Network::Undirected);
    net.populate(10000);
    net.erdos_renyi(5);

    cout << "days_ahead" << "," << "transmissibility1" << "," << "transmissibility2" << "," << "strain1" << "," << "strain2" << endl;

    for (int days_ahead = 0; days_ahead < 20; days_ahead ++){
        for (double transmissibility1 = 0.2; transmissibility1 <=  0.5; transmissibility1 += 0.05){
            for (double transmissibility2 = 0.2; transmissibility2 <= 0.5; transmissibility2 += 0.05){
                int i = 0;
                while ( i < 100){
                    // Choose and run simulation
                    Two_Strain_Percolation_Sim sim(&net);
                    sim.set_transmissibility1(transmissibility1);
                    sim.set_transmissibility2(transmissibility2);
                    sim.rand_infect(5, 0);
                    sim.run_headstart_simulation(days_ahead, 5);

                    if(sim.epidemic_size1() > 50){
                        cout << days_ahead << "," << transmissibility1 << "," << transmissibility2 << "," << sim.epidemic_size1() << "," << sim.epidemic_size2() << endl; 
                        i++;
                    }
                    
                    sim.reset();
                }
            }
        }
    }

    return 0;
}
