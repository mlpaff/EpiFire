#include "Deterministic_Network_TwoStrain_Sim.h"

int main() {

    double init_rho = log(0.000001);
    double max_rho = log(1.0);
    double step_size = (max_rho - init_rho) / 100;
    for(double rho1=init_rho; rho1 < max_rho; rho1=rho1+step_size){
        for(double rho2=init_rho; rho2 < max_rho; rho2=rho2+step_size){
            // Setup Initial y[] conditions for initialization
            double theta = 1.0;
            double pI1 = exp(rho1);
            double pI2 = exp(rho2);

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
                //vector<double> states = sim.get_state()
                //cout << exp(rho1) << " " << exp(rho2) << " ";
                //sim.printY();
                sim.step_simulation(1);
                i++;
            }
            cout << exp(rho1) << " " << exp(rho2) << " ";
            sim.printY();

        }
    }
    //sim.run_simulation();

    return 0;

}


