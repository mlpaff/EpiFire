#ifndef GIL_TWOSTRAIN_NET_SIM_H
#define GIL_TWOSTRAIN_NET_SIM_H

#include <stdlib.h>
#include <vector>
#include <iostream>
#include <queue>
#include "Utility.h"
#include "Network.h"

using namespace std;

class Event {
    public:
        double time;
        char type;
        Node* node;
        Event(const Event& o) {  time=o.time; type=o.type; node=o.node; }
        Event(double t, char e, Node* n) { time=t; type=e; node=n; }
        Event& operator=(const Event& o) { time=o.time; type=o.type; node=o.node;}
};

class compTime {
    public:
        bool operator() (const Event* lhs, const Event* rhs) const {
            return (lhs->time>rhs->time);
        }

        bool operator() (const Event& lhs, const Event& rhs) const {
            return (lhs.time>rhs.time);
        }
};


class Gillespie_TwoStrain_Network_Sim {
    public:
        // S -> Inf01, I02
        typedef enum {
            S, I01, I02, P1, P2, S12, S21, I12, I21, R, cumI1, cumI2, STATE_SIZE // STATE_SIZE must be last
        } stateType;
                                    // constructor
        Gillespie_TwoStrain_Network_Sim ( Network* net, double a1, double a2, double g1, double g2, double b1, double b2, double p1, double p2, int intro) {
            network = net;
            alpha1 = a1;
            alpha2 = a2;
            gamma1=g1;
            gamma2=g2;
            beta1 = b1;
            beta2=b2;
            phi1=p1;
            phi2=p2;
            intro_time = intro;
            reset();
        }

        Network* network;           // population
        double alpha1;
        double alpha2;
        double beta1;
        double beta2;
        double gamma1;
        double gamma2;
        double phi1;
        double phi2;
        int intro_time;

        // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        vector<int> state_counts;   
        double Now;                 // Current "time" in simulation

        MTRand mtrand;              // RNG

        void run_simulation(double duration) {
            double start_time = Now;
            int day = (int) Now;

            vector<Node*> nodes = rand_choose_nodes(1,2);
            add_event(intro_time, 'e',  nodes[0]);

            printStatus();
            while (next_event() and Now < start_time + duration) {
                if ((int) Now > day) {
                    //cout << fixed << setprecision(5);
                    printStatus();
                    day = (int) Now;
                }

                continue;
            }
        }

        void printStatus(){
            vector<int> vS1;
            vS1.push_back(S);
            vS1.push_back(S21);
            vector<int> vS2;
            vS2.push_back(S);
            vS2.push_back(S12);
            vector<int> vI1;
            vI1.push_back(I01);
            vI1.push_back(I21);
            vector<int> vI2;
            vI2.push_back(I02);
            vI2.push_back(I12);
            vector<double> meanSDEff = network->mean_eff_deg_states(vS2);

            cout << (int) Now << ", "
              << state_counts[S] << ", "
              << state_counts[I01] << ", " 
              << state_counts[I02] << ", " 
              << state_counts[P1] << ", " 
              << state_counts[P2] << ", " 
              << state_counts[S12] << ", " 
              << state_counts[S21] << ", " 
              << state_counts[I12] << ", " 
              << state_counts[I21] << ", " 
              << state_counts[R] << ", "
              << state_counts[cumI1] << ", "
              << state_counts[cumI2] <<  ", "
              << network->mean_deg_states(vS1) <<  ", "
              << network->mean_deg_states(vS2) <<  ", "
              << network->mean_deg_states(vI1) <<  ", "
              << network->mean_deg_states(vI2) << ", "
              << meanSDEff[0] << ", "
              << meanSDEff[1] << endl;
        }

        int current_epidemic_size1() {
            return state_counts[I01] + state_counts[I21];
        }

        int current_epidemic_size2() {
            return state_counts[I02] + state_counts[I12];
        }

        int reset() {
            Now = 0.0;
            vector<Node*> nodes = network->get_nodes();
            for (unsigned int i = 0; i < network->size(); i++) nodes[i]->set_state(S);

            state_counts.clear();
            state_counts.resize(STATE_SIZE, 0);
            state_counts[S] = network->size();
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        // choose n nodes without replacement
        vector<Node*> rand_choose_nodes (int n, int disease) {
            assert(n > -1 and n <= network->size());
            vector<Node*> nodes = network->get_nodes();
            vector<Node*> sample(n);
            vector<int> sample_ids(n);

            int all_susceptible = 0;
            while(all_susceptible==0){
                rand_nchoosek(network->size(), sample_ids, &mtrand);
                all_susceptible=1;
                for (unsigned int i = 0; i < sample_ids.size(); i++) {
                    if(disease==1){
                        if(nodes[ sample_ids[i] ]->get_state() != S &&
                                nodes[ sample_ids[i] ]->get_state() != S21){
                            all_susceptible = 0;
                        }
                    } else if (disease==2){
                        if(nodes[ sample_ids[i] ]->get_state() != S &&
                                nodes[ sample_ids[i] ]->get_state() != S12){
                            all_susceptible = 0;
                        }
                    }
                }
            }
            Node* node;
            for (unsigned int i = 0; i < sample_ids.size(); i++) {
                node = nodes[ sample_ids[i] ];
                sample[i] = node;
            }
            return sample;
        }


        void rand_infect(int n, int disease) {   // randomly infect k people
            vector<Node*> sample = rand_choose_nodes(n, disease);
            for (unsigned int i = 0; i < sample.size(); i++) {
                if(disease==1){
                    infect1(sample[i]);
                } else if(disease==2){
                    infect2(sample[i]);
                }
            }
            return;
        }

        void infect1(Node* node) {
            assert(state_counts[S] + state_counts[S12] + state_counts[S21] + state_counts[P2] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(I01);
                state_counts[S]--;
                state_counts[I01]++;
                current_state=I01;
            } else if(node->get_state() == S21){
                node->set_state(I21);
                state_counts[S21]--;  // decrement susceptible groupjj
                state_counts[I21]++;      // increment exposed group
                current_state=I21;
            } else if(node->get_state() == P2) {
                node->set_state(I21);
                state_counts[P2]--;  // decrement susceptible groupjj
                state_counts[I21]++;      // increment exposed group
                current_state=I21;
            }
            state_counts[cumI1]++;

            // time to recovery
            double Tr = rand_exp(gamma1, &mtrand) + Now;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double Tc = rand_exp(beta1, &mtrand) + Now;
                    if ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'c', neighbors[i]);   // potential transmission event
                    }
                }
            }

            if(current_state==I01){
                add_event(Tr, 'p', node);
                double Ts = Tr + rand_exp(alpha1, &mtrand);
                add_event(Ts, 't', node);
            } else if (current_state==I21){
                add_event(Tr, 'R', node);
            }
            return;
        }


        void infect2(Node* node) {
            assert(state_counts[S] + state_counts[S12] + state_counts[S21] + state_counts[P1] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(I02);
                state_counts[S]--;
                state_counts[I02]++;
                current_state=I02;
            } else if(node->get_state() == S12){
                node->set_state(I12);
                state_counts[S12]--;  // decrement susceptible groupjj
                state_counts[I12]++;      // increment exposed group
                current_state=I12;
            } else if(node->get_state() == P1) {
                node->set_state(I12);
                state_counts[P1]--;  // decrement susceptible groupjj
                state_counts[I12]++;      // increment exposed group
                current_state=I12;
            }
            state_counts[cumI2]++;

            // time to recovery
            double Tr = rand_exp(gamma2, &mtrand) + Now;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double Tc = rand_exp(beta2, &mtrand) + Now;
                    if ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'd', neighbors[i]);   // potential transmission event
                    }
                }
            }

            if(current_state==I02){
                add_event(Tr, 'q', node);
                double Ts = Tr + rand_exp(alpha2, &mtrand);
                add_event(Ts, 'u', node);
            } else if (current_state==I12){
                add_event(Tr, 'r', node);
            }
            return;
        }

        int next_event() {
            if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            Node* node = event.node;
            if (event.type == 'p') {
                if(node->get_state() == I01){
                    node->set_state(P1);
                    state_counts[I01]--;
                    state_counts[P1]++;
                }
            } else if (event.type == 'q') {
                if(node->get_state()==I02){
                    node->set_state(P2);
                    state_counts[I02]--;
                    state_counts[P2]++;
                }
            } else if (event.type == 't') {
                if(node->get_state()==P1){
                    node->set_state(S12);
                    state_counts[P1]--;
                    state_counts[S12]++;
                }
            } else if (event.type == 'u') {
                if(node->get_state()==P2){
                    node->set_state(S21);
                    state_counts[P2]--;
                    state_counts[S21]++;
                }
            } else if (event.type == 'r') {
                node->set_state(R);
                state_counts[I12]--;
                state_counts[R]++;
            } else if (event.type == 'R') {
                node->set_state(R);
                state_counts[I21]--;
                state_counts[R]++;
            } else if (event.type == 'c' ) {
                int node_state = node->get_state();
                if (node_state == S || node_state == S21){
                    infect1(node);
                } else if (node_state==P2 ){
                    if(mtrand.rand() < phi2){
                        infect1(node);
                    }
                }
            } else if (event.type == 'd' ) {
                int node_state = node->get_state();
                if (node_state == S || node_state == S12){
                    infect2(node);
                } else if (node_state==P1 ){
                    if(mtrand.rand() < phi1){
                        infect2(node);
                    }
                }
            } else if (event.type == 'e' ) {
                rand_infect(1, 2);
            } else {
                cerr << "Unknown event type encountered in simulator: " << event.type << "\nQuitting.\n";
            }
            return 1;
        }

        void add_event( double time, char type, Node* node) {
            EventQ.push( Event(time,type,node) );
            return;
        }

};
#endif
