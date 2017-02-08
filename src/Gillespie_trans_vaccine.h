#ifndef GIL_TRANS_VAX
#define GIL_TRANS_VAX

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
        //Event& operator=(const Event& o) { time=o.time; type=o.type; node=o.node;}
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


class Gillespie_trans_vaccine {
    public:
        typedef enum {
            S, I_vax, I_dis, R, cum_vax, cum_dis, STATE_SIZE // STATE_SIZE must be last
        } stateType;
                                    // constructor
        Gillespie_trans_vaccine ( Network* net, double g_vax, double g_dis, double b_vax, double b_dis, string vax_plan) {
            network = net;
            gamma_vax = g_vax;
            gamma_dis = g_dis;
            beta_vax = b_vax;
            beta_dis = b_dis;
            vax_type = vax_plan;
            reset();
        }

        Network* network;           // population
        double beta_vax;
        double beta_dis;
        double gamma_vax;
        double gamma_dis;
        string vax_type;

        // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        vector<int> state_counts;   
        double Now;                 // Current "time" in simulation

        MTRand mtrand;              // RNG

        void run_simulation(double duration, int init_vax, int init_infect) {
            
            if(vax_type == "trans"){
                rand_infect(init_vax, 1);
            } else if (vax_type == "rand"){

            } else {
                cerr << "Unknown vaccination type: " << vax_type << "\nQuitting.\n";
            }
            
            while (next_event()) { // Let all vaccination transmission occur
                //printStatus();
                continue;
            }
            
            Now = 0.0;
            double start_time = Now;
            int day = (int) Now;
            rand_infect(init_infect, 2);
            printStatus();
            while (next_event() and Now < start_time + duration) {
                if ((int) Now > day) {
                    printStatus();
                    day = (int) Now;
                }

                continue;
            }
        }

        void printStatus(){
            cout << (int) Now << ", "
              << state_counts[S] << ", "
              << state_counts[I_vax] << ", " 
              << state_counts[I_dis] << ", " 
              << state_counts[R] << ", "
              << state_counts[cum_vax] << ", "
              << state_counts[cum_dis] <<  endl;
        }

        int current_epidemic_size() {
            return state_counts[I_dis];
        }

        int current_vax_size() {
            return state_counts[I_vax];
        }

        void reset() {
            // Look this up
            Now = 0.0;
            vector<Node*> nodes = network->get_nodes();
            for (unsigned int i = 0; i < network->size(); i++) nodes[i]->set_state(S);

            state_counts.clear();
            state_counts.resize(STATE_SIZE, 0);
            state_counts[S] = network->size();
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
        }

        // choose n nodes without replacement
        vector<Node*> rand_choose_nodes (int n) {
            assert(n > -1 and n <= network->size());
            vector<Node*> nodes = network->get_nodes();
            vector<Node*> sample(n);
            vector<int> sample_ids(n);

            int all_susceptible = 0;
            while(all_susceptible==0){
                rand_nchoosek(network->size(), sample_ids, &mtrand);
                all_susceptible=1;
                for (unsigned int i = 0; i < sample_ids.size(); i++) {
                    if(nodes[ sample_ids[i] ]->get_state() != S){
                        all_susceptible = 0;
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
            vector<Node*> sample = rand_choose_nodes(n);
            for (unsigned int i = 0; i < sample.size(); i++) {
                if(disease==1){
                    infect_vax(sample[i]);
                } else if(disease==2){
                    infect_disease(sample[i]);
                }
            }
            return;
        }

        void infect_vax(Node* node) {
            assert(state_counts[S] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(I_vax);
                state_counts[S]--;
                state_counts[I_vax]++;
                state_counts[cum_vax]++;
            } else {
                cerr << "Vaccinating not a susceptible individual: " << node->get_state() << "\nQuitting.\n";
            }

            // time to recovery
            double Tr = rand_exp(gamma_vax, &mtrand) + Now;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double Tc = rand_exp(beta_vax, &mtrand) + Now;
                    while ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'v', neighbors[i]);   // potential transmission event
                        Tc = rand_exp(beta_vax, &mtrand) + Tc;
                    }
                }
            }
            add_event(Tr, 'r', node);

            return;
        }

        void infect_disease(Node* node) {
            assert(state_counts[S] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(I_dis);
                state_counts[S]--;
                state_counts[I_dis]++;
                state_counts[cum_dis]++;
            } else {
                cerr << "Infecting not a susceptible individual: " << node->get_state() << "\nQuitting.\n";
            }

            // time to recovery
            double Tr = rand_exp(gamma_dis, &mtrand) + Now;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double Tc = rand_exp(beta_dis, &mtrand) + Now;
                    while ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'i', neighbors[i]);   // potential transmission event
                        Tc = rand_exp(beta_dis, &mtrand) + Tc;
                    }
                }
            }
            add_event(Tr, 'R', node);

            return;
        }

        int next_event() {
            if ( EventQ.empty() ) return 0;
            Event event = EventQ.top(); // get the element
            EventQ.pop();               // remove from Q

            Now = event.time;           // advance time
            Node* node = event.node;
            
            if (event.type == 'r') {
                node->set_state(R);
                state_counts[I_vax]--;
                state_counts[R]++;
            } else if (event.type == 'R') {
                node->set_state(R);
                state_counts[I_dis]--;
                state_counts[R]++;
            } else if (event.type == 'i' ) {
                int node_state = node->get_state();
                if (node_state == S){
                    infect_disease(node);
                } 
            } else if (event.type == 'v' ) {
                int node_state = node->get_state();
                if (node_state == S){
                    infect_vax(node);
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
