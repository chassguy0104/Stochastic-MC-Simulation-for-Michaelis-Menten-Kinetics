/*Monte Carlo Simulation for Michaelis Menton Reaction -> Gillespie Algorithm




*/
#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <fstream>
using namespace std;

int main() {
    ofstream data("MCsimulation_data.txt");
    if (!data.is_open()) {
        cout << "Error opening file!" << endl;
        return -1;
    }
    data << "Trial\t\tTime\t\tS\t\tE\t\tES\t\tP" << endl;
    cout << endl;
    int total_trials = 100;
    long double max_time = 50.0;
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<long double> dist(0.0, 1.0);

    //Ensemble of Monte Carlo Trials
    for (int trial=1; trial<=total_trials; trial++) {

        double current_time = 0.0;
        int S=1000, E=100, ES=0, P=0;
        /*
         Here:
         S is the count of the substrate molecules
         E is the count of the enzyme molecules
         ES is the count of the enzyme-substrate complex molecules
         P is the count of the product molecules
         Michaelis-Menten kinetics relies on the Quasi-Steady-state assumption: S >> E
         S:E ratio being 10:1 provides a better visual plot than greater ratios (100:1 for instance)
         */
        double k1=0.002, k2=0.05, k3=0.1;
        /*
         Here:
         k1, k2, and k3 are the reaction rates of the 1st, 2nd and 3rd reactions
         As to why these particular values are chosen:

         */
        data <<  trial << "\t\t" << current_time << "\t\t" << S << "\t\t" << E << "\t\t" << ES << "\t\t" << P << endl;

        //Gillespie Algorithm : loops as long as time doesnt exceed max threshold AND molecular count of reactants is positive (reaction always proceeds forward)
        while (current_time < max_time && (S>0 || ES>0)) {
            //Propensities:
            double p1 = k1*S*E;
            double p2 = k2*ES;
            double p3 = k3*ES;
            double p_total = p1 + p2 + p3;

            if (p_total==0) break; //Nothing left to react

            double rand_1 = dist(gen);
            if (rand_1==0.0) {
                rand_1 = 1e-10; //to prevent log0 error
            }
            double rand_2 = dist(gen);

            double tau = (1.0/p_total)*log(1.0/rand_1); //time step -> time taken until next reaction starts
            current_time += tau;

            double threshold = rand_2 * p_total;
            if (threshold < p1) {
                //Binding reaction occurs: Reaction 1 is more inclined (higher propensity) to occur
                S--;
                E--;
                ES++;
            }
            else if (threshold < p1 + p2) {
                //Unbinding reaction occurs: Reaction 2 is more inclined to occur
                ES--;
                E++;
                S++;
            }
            else {
                //Catalysis (Reaction 3) is more inclined to occur
                ES--;
                P++;
                E++;
            }
            data << trial << "\t\t" << current_time << "\t\t" << S << "\t\t" << E << "\t\t" << ES << "\t\t" << P << endl;
        }

        cout << "Finished trial " << trial << " / " << total_trials << endl;


    }

    data.close();
    cout << "All data successfully saved to MCsimulation_data.txt" << endl;
    return 0;
}