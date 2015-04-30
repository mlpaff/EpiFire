#ifndef DETERMINISTIC_NETWORK_TWOSTRAIN_H
#define DETERMINISTIC_NETWORK_TWOSTRAIN_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_errno.h>
#include <iostream>
#include <math.h>
#include "DiffEq_Sim.h"

using namespace std;

class Deterministic_Network_TwoStrain_Sim : public DiffEq_Sim {

    private:
        const double ps_t0; //probability neighbor is initially susceptible
        const vector<double> deg_dist;
        const vector<double> deg_susc;
        const double beta1; // transmission rate 1
        const double beta2; // transmission rate 2
        const double gamma1; // recovery rate 1
        const double gamma2; // recovery rate 2


    public:
        Deterministic_Network_TwoStrain_Sim() : ps_t0(0.0), beta1(0.0),
            beta2(0.0), gamma1(0.0), gamma2(0.0){ nbins=7;}
        Deterministic_Network_TwoStrain_Sim(double ps_t0_param, double beta1_param, double beta2_param, double gamma1_param, double gamma2_param, vector<double> deg_dist_param, vector<double> deg_susc_param):
            ps_t0(ps_t0_param),
            beta1(beta1_param),
            beta2(beta2_param),
            gamma1(gamma1_param),
            gamma2(gamma2_param),
            deg_dist(deg_dist_param),
            deg_susc(deg_susc_param) {
                nbins=7;
            }
        ~Deterministic_Network_TwoStrain_Sim() {};

        void initialize( double theta, double pI1, double pI2) {
            y = new double[nbins];
            y[0] = theta;
            y[1] = pI1;
            y[2] = pI2;
            y[3] = pI1;
            y[4] = pI2;
            y[5] = 0.0;
            y[6] = 0.0;
        }

        double current_susceptible() { return g( y[0] ); }
        double current_infectious_1() { return y[3]; }
        double current_infectious_2() { return y[4]; }
        double current_recovered_1() { return y[5]; }
        double current_recovered_2() { return y[6]; }

        double current_all_infectious() {
            return current_infectious_1() + current_infectious_2();
        }

        double current_all_recovered() {
            return current_recovered_1() + current_recovered_2();
        }


        double g( double theta) {
            double val = 0;
            for (unsigned int i = 0; i<deg_dist.size(); i++) {
                val += deg_dist[i] * pow(theta,i);
            }
            return val;
        }

        double dg( double theta) {
            double val = 0;
            for (unsigned int i = 1; i<deg_dist.size(); i++) {
                val += (i) * deg_dist[i] * pow(theta,i-1);
            }
            return val;
        }

        double ddg( double theta) {
            double val = 0;
            for (unsigned int i = 2; i<deg_dist.size(); i++) {
                val += (i) * (i-1) *deg_dist[i] * pow(theta,i-2);
            }
            return val;
        }


        void derivative(double const y[], double dydt[]) {
            const double theta = y[0];
            const double pI1 = y[1];
            const double pI2 = y[2];
            const double I1 = y[3];
            const double I2 = y[4];
            const double R1 = y[5];
            const double R2 = y[6];

            dydt[0] = - beta1*pI1 - beta2*pI2;                  // dtheta.dt
            dydt[1] = -(beta1+gamma1)*pI1 + beta1*pI1*ps_t0*ddg(theta)/dg(1);      // dpI1.dt
            dydt[2] = -(beta2+gamma2)*pI2 + beta2*pI2*ps_t0*ddg(theta)/dg(1);      // dpI2.dt
            dydt[3] = beta1*pI1*dg(theta) - gamma1*I1; // dI1.dt
            dydt[4] = beta2*pI2*dg(theta) - gamma2*I2; // dI2.dt
            dydt[5] = gamma1*I1;   //dR.dt
            dydt[6] = gamma2*I2;   //dR.dt
        }

};

#endif

