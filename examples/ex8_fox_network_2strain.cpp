#include "Deterministic_Network_TwoStrain_Sim.h"

int main() {


    // Setup Initial y[] conditions for initialization
    double pop_size = 10000.0;
    double theta = 1.0;
    double pI1 = 1/pop_size;
    double pI2 = 1/pop_size;

    // Setup model parameters
    double ps_t0 = 1 - pI1 - pI2;
    double beta1 = 1.2;
    double beta2 = 0.2;
    double gamma1 = 4.0;
    double gamma2 = 0.25;

    // setup uniform network degree 1-9 as
    // joel Miller phys Rev paper 2013 Fig 3
    vector<double> degree_dist;
    degree_dist.push_back(0.0);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);
    degree_dist.push_back(0.11);

    suscs

    Deterministic_Network_TwoStrain_Sim sim(ps_t0, beta1, beta2, gamma1, gamma2, degree_dist, degree_dist);
    sim.initialize(theta, pI1, pI2);
    int i = 0;
    while (i < 100) {
        sim.printY();
        sim.step_simulation(1);
        i++;
    }
    sim.printY();

    //sim.run_simulation();

    return 0;
}


