#include "Gill_SEIR_TwoStrain_Seasonal_Network.h"
// runs a seasonally forced gillespie simulation on a exponential network
int main(int argc,char *argv[]) {
    if(argc != 9) {
            printf("Wrong Number of Arguments\n");
            exit(0);
    }

    double alpha1, alpha2, beta1_max, beta2_max, gamma1, gamma2, phi1, phi2, eta1, eta2;
    alpha1 = alpha2 = atof( argv[1] );
    beta1_max = atof( argv[2] );  
    beta2_max = atof( argv[3] );
    eta1 = eta2 = 1.0/2.62;
    gamma1 = gamma2 = 1.0/3.38;
    phi1 = phi2 = atof( argv[4]);
    int intro_time, start_ind, shift;
    intro_time = atoi( argv[5] );
    start_ind = atoi( argv[6] );
    shift = atoi( argv[7] );
    string network_type = argv[8];
    
    int num_reps = 5000;

    Network net = Network("gillespie toy", Network::Undirected);
    net.populate(10000);
    for(int i =1; i <= num_reps; i++){
        net.clear_edges();
        if(network_type == "exp"){
            net.rand_connect_exponential(0.06453487);
        } else if(network_type == "unif"){
            vector<int> degrees(10000, 16);
            net.rand_connect_explicit(degrees);
        } else{
            cerr << "Unrecognized network type" << endl;
        }
        // cout << net.mean_deg() << endl;
        Gillespie_SEIR_TwoStrain_Network sim(&net, alpha1, alpha2, eta1, eta2, gamma1, gamma2, beta1_max, beta2_max, phi1, phi2, intro_time, start_ind, shift);
        cout << "Simulation number: " << i << endl;
        sim.reset();
        sim.rand_infect(5, 1);
        sim.run_simulation(10000.0);
    }
    return 0;
}
