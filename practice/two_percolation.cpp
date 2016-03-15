#include <Two_Strain_Percolation_Sim.h>

int main() {
    
        // Construct Network
        Network net("name", Network::Undirected);
        net.populate(10000);
        net.erdos_renyi(5);

    cout << "transmissibility" << ","<< "expected_r0" << ","  << "epidemic_size" << endl;
    for(double t = 0; t<1; t += 0.05){
        vector<int> epi_sizes;
        double expected_r0;
        for (int i = 0; i < 1000; i++){
            // Choose and run simulation
            Two_Strain_Percolation_Sim sim(&net);
            sim.set_transmissibility1(t);
            sim.rand_infect(5, 0);
            expected_r0 = sim.expected_R0_1();
            sim.run_simulation();
            epi_sizes.push_back(sim.epidemic_size1());
            // cout << sim.epidemic_size1() << "\t" << sim.epidemic_size2() << endl;
            sim.reset();
        }   
        cout << t << "," << expected_r0 << ","  << mean(epi_sizes) << endl;
    }

    return 0;
}
