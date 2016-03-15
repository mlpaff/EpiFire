#include "Deterministic_Network_TwoStrain_Sim.h"

int main() {

    // Setup Initial y[] conditions for initialization
    double theta = 1.0;
    double pI1 = 0.001;
    double pI2 = 0.001;

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

    vector<double> degree_susceptibility;
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);
    degree_susceptibility.push_back(ps_t0);


    Deterministic_Network_TwoStrain_Sim sim(ps_t0, beta1, beta2, gamma1, gamma2, degree_dist, degree_susceptibility);
    sim.initialize(theta, pI1, pI2);
    int i = 0;
    while (i < 99) {
        sim.printY();
        sim.step_simulation(1);
        i++;
    }
    sim.printY();
    cout << sim.growth_rate_1() << " " << sim.growth_rate_2() << " " << sim.ddg(1.0) << " " << sim.dg(1.0) <<endl;
    return 0;
}


