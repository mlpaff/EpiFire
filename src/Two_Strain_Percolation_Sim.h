#ifndef TWO_STRAIN_PERCOL_SIMULATOR_H
#define TWO_STRAIN_PERCOL_SIMULATOR_H

#include "Simulator.h"

class Two_Strain_Percolation_Sim: public Simulator
{
    protected:
        vector<Node*> infected1;
        vector<Node*> infected2;
        vector<Node*> recovered1;
        vector<Node*> recovered2;

    public:
        typedef enum {           //Whatever is equal to zero is the default state
            S=0, I1=1, I2 = 2, R1=-1, R2=-2
        } stateType;
        float T1, T2;                 // transmissibiltiy, e.g. P{infection spreading along a given edge that connects an infected and a susceptible}

        Two_Strain_Percolation_Sim():Simulator() {T1=0; T2=0;};
        Two_Strain_Percolation_Sim(Network* net):Simulator(net) {T1=0; T2=0;};
        ~Two_Strain_Percolation_Sim() {};

        void set_transmissibility1(double t) { this->T1 = t; }
        void set_transmissibility2(double t) { this->T2 = t; }


        double expected_R0_1 () {
            //assert(T != NULL);
            return T1 / calc_critical_transmissibility();
        }

        double expected_R0_2 () {
            //assert(T != NULL);
            return T2 / calc_critical_transmissibility();
        }

        vector<Node*> rand_expose (int n, int m) {
            int k = n+m;
            assert(k >= 0);
            if (k == 0) return vector<Node*>(0);
            vector<Node*> patients_zero = rand_choose_nodes(k);
            stateType state = I1;
            vector<Node*>* tmp_infected = &infected1;
            for (int i = 0; i < k; i++) {
                if (i>=n) {
                    state = I2;
                    tmp_infected = &infected2;
                }
                if (patients_zero[i]->get_state() == S) {
                    patients_zero[i] -> set_state( state );
                    (*tmp_infected).push_back(patients_zero[i]);
                }
            }
                
            return patients_zero;
        }

        vector<Node*> rand_infect (int n, int m) {
            assert((n + m) > 0);
            vector<Node*> patients_zero = rand_choose_nodes(n+m);
            for (int i = 0; i < n; i++) {
                patients_zero[i] -> set_state( I1 );
                infected1.push_back(patients_zero[i]);
            };       
            for (int i = 0; i < m; i++) {
                patients_zero[i + n] -> set_state( I2 );
                infected2.push_back(patients_zero[i + n]);
            };
            return patients_zero;
        }

        void step_simulation1 () {
            assert(infected1.size() > 0);
            time++;
            //cerr << "\t" << infected.size() << endl;
            vector<Node*> new_infected_1;
            for (unsigned int i = 0; i < infected1.size(); i++) {
                Node* inode = infected1[i];
                vector<Node*> neighbors = inode->get_rand_neighbors();

                for (unsigned int j = 0; j < neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( test->get_state() == S && mtrand->rand() < T1 ) {
                        test->set_state( I1 );
                        new_infected_1.push_back( test );
                    }
                }
                inode->set_state( R1 );
                recovered1.push_back( inode );
            }
            
            infected1 = new_infected_1;
        }

        void step_simulation2 () {
            assert(infected2.size() > 0);
            time++;
            //cerr << "\t" << infected.size() << endl;
            vector<Node*> new_infected_2;
            for (unsigned int i = 0; i < infected2.size(); i++) {
                Node* inode = infected2[i];
                vector<Node*> neighbors = inode->get_rand_neighbors();

                for (unsigned int j = 0; j < neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( test->get_state() == S && mtrand->rand() < T1 ) {
                        test->set_state( I2 );
                        new_infected_2.push_back( test );
                    }
                }
                inode->set_state( R2 );
                recovered2.push_back( inode );
            }
            
            infected2 = new_infected_2;
        }

        void step_simulation () {
            assert(infected1.size() > 0 || infected2.size() > 0);
            time++;
            //cerr << "\t" << infected.size() << endl;
            vector<Node*> new_infected_1;
            vector<Node*> new_infected_2;
            for (unsigned int i = 0; i < infected1.size(); i++) {
                Node* inode = infected1[i];
                vector<Node*> neighbors = inode->get_rand_neighbors();

                for (unsigned int j = 0; j < neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( test->get_state() == S && mtrand->rand() < T1 ) {
                        test->set_state( I1 );
                        new_infected_1.push_back( test );
                    }
                }
                inode->set_state( R1 );
                recovered1.push_back( inode );
            }
            

            for (Node* inode: infected2) {
                vector<Node*> neighbors = inode->get_rand_neighbors();

                for (unsigned int j = 0; j < neighbors.size(); j++) {
                    Node* test = neighbors[j];
                    if ( mtrand->rand() < T2 ) {
                        if(test->get_state() == S){
                            test->set_state( I2 );
                            new_infected_2.push_back( test );
                        }
                        else if(test->get_state() == I1 && mtrand->rand() < (T2 / (T1+T2)) ){
                            test->set_state( I2 );
                            new_infected_2.push_back( test );
                            new_infected_1.erase(std::find(new_infected_1.begin(), new_infected_1.end(), test) );
                        }
                    }
                }
                inode->set_state( R2 );
                recovered2.push_back( inode );
            }
            
            infected1 = new_infected_1;
            infected2 = new_infected_2;
        }

        void run_simulation() {
            while (infected1.size() > 0 || infected2.size() > 0) {
                step_simulation();
            }
        }

        void run_headstart_simulation(int days, int num_expose2) {
            for (int day = 0; day < days; day++) {
                if(infected1.size() > 0) step_simulation1(); 
            }
            rand_expose(0, num_expose2);

            while (infected1.size() > 0 || infected2.size() > 0) {
                step_simulation();
            }
        }

        int count_infected1() {
            return infected1.size();
        }
        int count_infected2() {
            return infected2.size();
        }

        int epidemic_size() { return recovered1.size() + recovered2.size(); }
        int epidemic_size1() { return recovered1.size(); }
        int epidemic_size2() { return recovered2.size(); }

        void reset() {
            reset_time();

            set_these_nodes_to_state(infected1, S);
            infected1.clear();

            set_these_nodes_to_state(recovered1, S);
            recovered1.clear();

            set_these_nodes_to_state(infected2, S);
            infected2.clear();

            set_these_nodes_to_state(recovered2, S);
            recovered2.clear();
        }

        void summary() {
            cerr << "Network size: " << net->size();
            cerr << "\tTransmissibility: " << T1 << "\t" << T2;
            cerr << "\tEpidemic size: " << recovered1.size() << "\t" << recovered2.size() << "\n\n";

        }
};
#endif
