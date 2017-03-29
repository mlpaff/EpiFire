#ifndef GIL_SEIR_TWOSTRAIN_SEAS_NET_H
#define GIL_SEIR_TWOSTRAIN_SEAS_NET_H

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


class Gillespie_SEIR_TwoStrain_Network {
    public:
        // S -> Inf01, I02
        typedef enum {
            S, E01, E02, I01, I02, P1, P2, S12, S21, E12, E21, I12, I21, R, cumI1, cumI2, STATE_SIZE // STATE_SIZE must be last
        } stateType;
                                    // constructor
        Gillespie_SEIR_TwoStrain_Network ( Network* net, double a1, double a2, double e1, double e2, 
                                            double g1, double g2, double b1, double b2, double p1, double p2, int intro,
                                            int strt_ind, int shft) {
            network = net;
            alpha1 = a1;
            alpha2 = a2;
            eta1 = e1;
            eta2=e2;
            gamma1=g1;
            gamma2=g2;
            beta1_max=b1;
            beta2_max=b2;
            phi1=p1;
            phi2=p2;
            intro_time = intro;
            seas_start_ind = strt_ind;
            shift = shft;
            reset();
        }

        Network* network;           // population
        double alpha1;
        double alpha2;
        double beta1_max;
        double beta2_max;
        double eta1;
        double eta2;
        double gamma1;
        double gamma2;
        double phi1;
        double phi2;
        int intro_time;
        int seas_start_ind;
        int shift;
        double psi; // network calculation



        // event queue
        priority_queue<Event, vector<Event>, compTime > EventQ;
        vector<int> state_counts;   
        double Now;                 // Current "time" in simulation
        vector<double> hum_data;

        MTRand mtrand;              // RNG

        void run_simulation(double duration) {
            double start_time = Now;
            int day = (int) Now;

            vector<Node*> nodes = rand_choose_nodes(1,2);
            add_event(intro_time, 'e',  nodes[0]);

            printStatus();
            //printNetStatus();
            while (next_event() and Now < start_time + duration) {
                if ((int) Now > day) {
                    printStatus();
                    //printNetStatus();
                    day = (int) Now;
                }

                continue;
            }
        }

        void printNetStatus(){
            vector<int> vS2;
            vS2.push_back(S);
            vS2.push_back(S12);

            vector< vector<int> > stateDeg = network->get_states_by_degree();
            cout << (int) Now << ", ";
            for(int i=0; i<stateDeg.size(); i++){
                int suscCount =0;
                for(int j=0; j<stateDeg[i].size(); j++){
                    for(int k=0; k < vS2.size();k++){
                        if(stateDeg[i][j]==vS2[k]){
                            suscCount ++;
                        }
                    }
                }
                if(i==(stateDeg.size()-1)) {
                    cout << suscCount << endl;
                }else{
                    cout << suscCount << ", ";
                }
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
              << state_counts[E01] << ", " 
              << state_counts[E02] << ", " 
              << state_counts[I01] << ", " 
              << state_counts[I02] << ", " 
              << state_counts[P1] << ", " 
              << state_counts[P2] << ", " 
              << state_counts[S12] << ", " 
              << state_counts[S21] << ", " 
              << state_counts[E12] << ", " 
              << state_counts[E21] << ", " 
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
            // Look this up
            Now = 0.0;
            psi = ( (network->mean_sq_deg() - network->mean_deg()) / network->mean_deg() );

            vector<Node*> nodes = network->get_nodes();
            for (unsigned int i = 0; i < network->size(); i++) nodes[i]->set_state(S);

            state_counts.clear();
            state_counts.resize(STATE_SIZE, 0);
            state_counts[S] = network->size();
            EventQ = priority_queue<Event, vector<Event>, compTime > ();
            double hum_data_arr[365] = {0.003198947,0.00322633,0.003289128,0.003325971,0.003298366,
                0.00331145,0.003292259,0.003364593,0.003393615,0.003368119,0.003315297,0.003290407,0.003194007,0.003025532,
                0.002958344,0.002994219,0.003076709,0.003119807,0.003087744,0.002999556,0.002998688,0.003030747,0.003130397,
                0.003183583,0.003165547,0.00315221,0.00314404,0.003147307,0.00319955,0.003217159,0.003121045,0.003141987,
                0.003169503,0.003193127,0.00318428,0.003191315,0.003147461,0.003159119,0.003133069,0.003093632,0.003064768,
                0.003051333,0.003120598,0.003208409,0.003247492,0.003291129,0.003300963,0.003227149,0.003255415,0.003289257,
                0.003357683,0.003439811,0.003462622,0.003387482,0.003329617,0.00327451,0.003238261,0.003271157,0.003295243,
                0.003397928,0.003485081,0.003426192,0.003400088,0.003479221,0.003540497,0.003653592,0.003762976,0.003889279,
                0.00387103,0.00391348,0.004016564,0.004086288,0.004089897,0.004137052,0.00414938,0.004075651,0.00400936,
                0.004093356,0.004180954,0.004141574,0.004092302,0.004059425,0.004084505,0.004166909,0.004230767,0.004269329,
                0.004401611,0.004503554,0.004534535,0.004613848,0.004632451,0.004605263,0.004543125,0.004463418,0.00446199,
                0.004486839,0.004615195,0.004753607,0.004767002,0.004800323,0.004896122,0.004942022,0.005004858,0.004971983,
                0.004990521,0.005038833,0.005127721,0.005215592,0.005300411,0.005360101,0.005446078,0.005434538,0.0053848,
                0.005418233,0.005414044,0.005474047,0.00554793,0.00553087,0.005513124,0.005630799,0.005789668,0.005929235,
                0.005937052,0.005960824,0.006081983,0.006273586,0.006432604,0.006552961,0.00654488,0.00649442,0.006450798,
                0.006428501,0.006418515,0.006548398,0.006645633,0.006760629,0.006786966,0.006872082,0.007043075,0.007052734,
                0.007087889,0.007205641,0.007371201,0.007475906,0.007637161,0.007780515,0.007952893,0.008095483,0.00808607,
                0.008028129,0.008119466,0.008186571,0.008203834,0.00826218,0.008286997,0.008438398,0.008470981,0.008421386,
                0.008516449,0.008578557,0.008678719,0.008886438,0.009002821,0.009026997,0.009033907,0.009034974,0.009140992,
                0.009167129,0.009150722,0.009178349,0.009227477,0.009338501,0.009506295,0.009575255,0.009530726,0.009657782,
                0.009846878,0.009947677,0.009948041,0.009938318,0.009892233,0.009847438,0.009929744,0.010016596,0.01015316,
                0.01021659,0.01023933,0.010283983,0.010340935,0.010414548,0.010526102,0.010624379,0.010705534,0.010733239,
                0.010688241,0.010631048,0.010689278,0.010768013,0.010854637,0.010847611,0.010810207,0.010871297,0.010890763,
                0.010864803,0.010832833,0.01086835,0.010946145,0.011010552,0.011040965,0.011011997,0.010968744,0.01102418,
                0.011090011,0.010980811,0.010887562,0.010866904,0.010813023,0.010746208,0.010728705,0.010796219,0.010830876,
                0.010782498,0.010663285,0.010636888,0.010569036,0.010496687,0.010483796,0.010436851,0.010346209,0.010364968,
                0.010416951,0.010354168,0.01030635,0.010257482,0.010125777,0.010052367,0.010093394,0.010192663,0.01018844,
                0.01022242,0.010220901,0.010105029,0.009966109,0.009896147,0.009817427,0.009633393,0.009451563,0.009395076,
                0.009428996,0.009352573,0.009326428,0.009245499,0.009123207,0.008934504,0.008782848,0.008663736,0.00856986,
                0.008528717,0.008422733,0.008408076,0.008403253,0.008403532,0.008332793,0.008219307,0.008124004,0.008029764,
                0.007968998,0.007934937,0.007850213,0.007742778,0.00762558,0.007394711,0.0073015,0.007295495,0.007226354,
                0.00712183,0.006989496,0.006869294,0.00687845,0.006827904,0.006713933,0.006716971,0.006739723,0.006650732,
                0.006513685,0.006432885,0.006363315,0.006195925,0.006062685,0.006075981,0.006084838,0.006118457,0.006013904,
                0.005943209,0.005883896,0.005794794,0.005760793,0.005690344,0.005614966,0.005481069,0.005409325,0.005347253,
                0.005310586,0.005339769,0.005280994,0.005170782,0.005088169,0.005087738,0.005155757,0.005111778,0.004931265,
                0.004835775,0.004860791,0.004837885,0.004698099,0.004593282,0.004595591,0.004631508,
                0.004626624,0.004603629,0.004520982,0.004355186,0.004273496,0.004198169,0.004178147,0.004133834,
                0.004028545,0.003975883,0.003948585,0.004000508,0.003902346,0.003827881,0.003877692,0.00388316,
                0.003805025,0.003708332,0.003664678,0.00363753,0.00360868,0.003494694,0.003459791,0.003443582,
                0.003504107,0.003560239,0.003605335,0.003597749,0.00366062,0.003603812,0.003514914,0.003407986,
                0.003437623,0.003418602,0.003331193,0.003400198,0.003463376,0.003531663,0.003517699,0.003463488,
                0.003409076,0.00333295,0.003302318,0.003286818,0.003288187,0.003244882,0.003184149};
            for(int i =0; i < 365; i++){
                hum_data.push_back(hum_data_arr[i]);
            }
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


        double get_seasonal_beta(double time_to_add_to_now, double beta_max, double gamma){
            double rnot_max = beta_max / (beta_max + gamma) * psi;
            double timept = Now + time_to_add_to_now + seas_start_ind - shift;
            if(timept <0){
                timept = 365 + timept;
            }
            int hum_index =  (int) (fmod(timept, 365.0) - 1.0);
            if(hum_index == -1 ){
                hum_index = 364;
            }

            double rnot_t = exp(-180.0 * hum_data[hum_index] + log(rnot_max - 0.1)) + 0.1;
            //cout << Now << " " << time_to_add_to_now << " " << hum_index << " " << hum_data[hum_index] << " " << rnot_t << endl;
            return((rnot_t * gamma / psi) / (1 - rnot_t / psi));
        }

        void infect1(Node* node) {
            assert(state_counts[S] + state_counts[S12] + state_counts[S21] + state_counts[P2] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(E01);
                state_counts[S]--;
                state_counts[E01]++;
                current_state=E01;
            } else if(node->get_state() == S21){
                node->set_state(E21);
                state_counts[S21]--;  // decrement susceptible groupjj
                state_counts[E21]++;      // increment exposed group
                current_state=E21;
            } else if(node->get_state() == P2) {
                node->set_state(E21);
                state_counts[P2]--;  // decrement susceptible groupjj
                state_counts[E21]++;      // increment exposed group
                current_state=E21;
            } else {
                cerr << "Unknown state type encountered in infection: " << node->get_state() << "\nQuitting.\n";
            }
           // state_counts[cumI1]++;

            // time to infectious
            double Ti = rand_exp(eta1, &mtrand) + Now;
            if(current_state==E01){
                add_event(Ti, 'f', node);
            } else if(current_state==E21){
                add_event(Ti, 'g', node);
            }

            // time to recovery
            double Tr = rand_exp(gamma1, &mtrand) + Ti;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double beta1 = get_seasonal_beta(Ti, beta1_max, gamma1);
                    double Tc = rand_exp(beta1, &mtrand) + Ti;
                    while ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'c', neighbors[i]);   // potential transmission event
                        Tc = rand_exp(beta1, &mtrand) + Tc;
                    }
                }
            }

            if(current_state==E01){
                add_event(Tr, 'p', node);
                double Ts = Tr + rand_exp(alpha1, &mtrand);
                add_event(Ts, 't', node);
            } else if (current_state==E21){
                add_event(Tr, 'R', node);
            }
            return;
        }


        void infect2(Node* node) {
            assert(state_counts[S] + state_counts[S12] + state_counts[S21] + state_counts[P1] > 0);
            int current_state;
            if(node->get_state() == S){
                node->set_state(E02);
                state_counts[S]--;
                state_counts[E02]++;
                current_state=E02;
            } else if(node->get_state() == S12){
                node->set_state(E12);
                state_counts[S12]--;  // decrement susceptible groupjj
                state_counts[E12]++;      // increment exposed group
                current_state=E12;
            } else if(node->get_state() == P1) {
                node->set_state(E12);
                state_counts[P1]--;  // decrement susceptible groupjj
                state_counts[E12]++;      // increment exposed group
                current_state=E12;
            } else {
                cerr << "Unknown state type encountered in infection 2: " << node->get_state() << "\nQuitting.\n";
            }
            //state_counts[cumI2]++;

            // time to infectious
            double Ti = rand_exp(eta2, &mtrand) + Now;
            if(current_state==E02){
                add_event(Ti, 'h', node);
            } else if(current_state==E12){
                add_event(Ti, 'i', node);
            }

            // time to recovery
            double Tr = rand_exp(gamma2, &mtrand) + Ti;
            // time to next contact
            vector<Node*> neighbors = node->get_neighbors();
            if(neighbors.size() > 0){
                for(int i=0; i<neighbors.size();i++){
                    double beta2 = get_seasonal_beta(Ti, beta2_max, gamma2);
                    double Tc = rand_exp(beta2, &mtrand) + Ti;
                    while ( Tc < Tr ) {     // does contact occur before recovery?
                        add_event(Tc, 'd', neighbors[i]);   // potential transmission event
                        Tc = Tc + rand_exp(beta2, &mtrand);
                    }
                }
            }

            if(current_state==E02){
                add_event(Tr, 'q', node);
                double Ts = Tr + rand_exp(alpha2, &mtrand);
                add_event(Ts, 'u', node);
            } else if (current_state==E12){
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
            } else if (event.type == 'f') {
                if(node->get_state()==E01){
                    node->set_state(I01);
                    state_counts[I01]++;
                    state_counts[E01]--;
                    state_counts[cumI1]++;
                }
            } else if (event.type == 'g') {
                if(node->get_state()==E21){
                    node->set_state(I21);
                    state_counts[E21]--;
                    state_counts[I21]++;
                    state_counts[cumI1]++;
                }
            } else if (event.type == 'h') {
                if(node->get_state()==E02){
                    node->set_state(I02);
                    state_counts[E02]--;
                    state_counts[I02]++;
                    state_counts[cumI2]++;
                }
            } else if (event.type == 'i') {
                if(node->get_state()==E12){
                    node->set_state(I12);
                    state_counts[E12]--;
                    state_counts[I12]++;
                    state_counts[cumI2]++;
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
