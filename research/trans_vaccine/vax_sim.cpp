#include "Gillespie_trans_vaccine.h"

int main(int argc,char *argv[]) {
    if(argc != 10) {
            printf("Wrong Number of Arguments\n");
            exit(0);
    }

    //double alpha1 = argv[1]-0.0;
    double beta_vax, beta_dis, gamma_vax, gamma_dis;
    int init_inf, init_vax;

    beta_vax = atof( argv[1] ); // vaccine transmissibility
    beta_dis = atof( argv[2] ); // disease transmissibility
    gamma_vax = atof( argv[3] ); // vaccine recovery rate
    gamma_dis = atof( argv[4] ); // disease recovery rate
    init_vax = atoi( argv[5] ); // number initially vaccinated
    init_inf = atoi( argv[6] ); // number initially infected

    string network_type = argv[7] ;
    string vax_type = argv[8];
    int num_reps = atoi( argv[9]);

    Network net = Network("gillespie toy", Network::Undirected);
    net.populate(10000);
    for(int i =1; i <= num_reps; i++){
        net.clear_edges();
        if(network_type == "er"){
            net.erdos_renyi(10);
        } else if(network_type == "exp"){
            net.rand_connect_exponential(0.1);
        } else if(network_type == "unif"){
            vector<int> degrees(10000, 10);
            net.rand_connect_explicit(degrees);
        } else{
            cerr << "Unrecognized network type" << endl;
        }
        Gillespie_trans_vaccine sim(&net, gamma_vax, gamma_dis, beta_vax, beta_dis, vax_type);
        cout << "Simulation number: " << i << endl;
        sim.reset();
        sim.run_simulation(10000.0, init_vax, init_inf);
    }
    return 0;
}


