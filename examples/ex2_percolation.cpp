#include <Percolation_Sim.h>

int main() {
    
        // Construct Network
        Network net("name", false);
        net.populate(10000);
        net.erdos_renyi(5);

    for (int i = 0; i < 10000; i++){
        // Choose and run simulation
        Percolation_Sim sim(&net);
        sim.set_transmissibility(0.2);
        sim.rand_infect(10);
        sim.run_simulation();
        cout << sim.epidemic_size() << endl;
        sim.reset();
    }
    return 0;
}
