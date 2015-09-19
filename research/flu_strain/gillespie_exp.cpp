#include "Gillespie_TwoStrain_Network_Sim.h"

int main(int argc,char *argv[]) {
    if(argc != 5) {
            printf("Wrong Number of Arguments\n");
            exit(0);
    }

    //double alpha1 = argv[1]-0.0;
    double alpha1, alpha2, beta1, beta2, gamma1, gamma2, phi1, phi2;
    alpha1 = alpha2 = atof( argv[1] );
    //beta1 = 0.003383833;
    beta1 = 0.043;
    beta2 = atof( argv[2] );
    gamma1 = gamma2 = 1.0/5.0;
    phi1 = phi2 = atof( argv[4]);
    int intro_time;
    intro_time = atoi( argv[3] );

    int num_reps = 1000;

    Network net = Network("gillespie toy", Network::Undirected);
    net.populate(10000);
    for(int i =1; i <= num_reps; i++){
        net.clear_edges();
        //net.rand_connect_powerlaw(2.0, 689.0);
        net.rand_connect_exponential(0.222);
        //cout << net.mean_deg() << endl;
        Gillespie_TwoStrain_Network_Sim sim(&net, alpha1, alpha2, gamma1, gamma2, beta1, beta2, phi1, phi2, intro_time);
        cout << "Simulation number: " << i << endl;
        sim.reset();
        sim.rand_infect(10, 1);
        sim.run_simulation(10000.0);
    }
    return 0;
}


