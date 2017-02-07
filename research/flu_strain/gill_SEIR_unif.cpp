#include "Gillespie_SEIR_TwoStrain_Network.h"
// runs a gillespie simulation on a uniform network
int main(int argc,char *argv[]) {
    if(argc != 5) {
            printf("Wrong Number of Arguments\n");
            exit(0);
    }

    double alpha1, alpha2, beta1, beta2, gamma1, gamma2, phi1, phi2, eta1, eta2;
    alpha1 = alpha2 = atof( argv[1] );
    beta1 = 0.03024476;    //2008 season fit
    //beta1 = 0.04063644;  //2003 season fit
    beta2 = atof( argv[2] );
    eta1 = eta2 = 1.0/2.62;
    gamma1 = gamma2 = 1.0/3.38;
    phi1 = phi2 = atof( argv[4]);
    int intro_time;
    intro_time = atoi( argv[3] );

    int num_reps = 1;

    Network net = Network("gillespie toy", Network::Undirected);
    net.populate(10000);
    vector<int> degrees(10000, 16);
    //cout << net.get_biggest_component().size() << endl;
    
    for(int i =1; i <= num_reps; i++){
        net.clear_edges();
        net.rand_connect_explicit(degrees);
        // cout << net.mean_deg()<<endl;
        Gillespie_SEIR_TwoStrain_Network sim(&net, alpha1, alpha2, eta1, eta2, gamma1, gamma2, beta1, beta2, phi1, phi2, intro_time);
        cout << "Simulation number: " << i << endl;
        sim.reset();
        sim.rand_infect(5, 1);
        sim.run_simulation(10000.0);
    }
    return 0;
}


