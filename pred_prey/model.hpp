#ifndef PRED_PREY_MODEL_HPP_INCLUDED 
#define PRED_PREY_MODEL_HPP_INCLUDED

#include <random>
#include <torus/torus.hpp>
#include <torus/torus_grid.hpp>
#include <torus/torus_hrtree.hpp>
#include <fstream>
#include <string>
#include <iostream>
#include "ann2.hpp"

namespace model {
    //*** parameters ************************************************************
    using namespace ann;
    using default_activation = activation::rtlu;
    using Brain = Network<double,
        Layer<Neuron<2, default_activation>, 1>
    >;
    //using Brain = Network<double,
    //    Layer<Neuron<3, default_activation>, 3>, 
    //    Layer<Neuron<3, default_activation>, 1>
    //>;


    struct Param
    {
        // simulation parameters
        const int generationTime = 50;
        const int maxTime = 300'000;// 200'000;
        size_t GS = 500;         // grid size
        const float randomWalk = 2.0f;

        // initial values
        size_t N_prey = 15'000;   // 6'000  
        size_t N_pred = 1000; // 400   
        const double resourceDensity = 0.05;
        const double initResourcesDensity = 0.05;

        // resources parameters
        const double grow_back = 0.005;  // resource grow back (increment)
        
        // prey parameters
        const float prey_sr = 4.0f; // search radius [expressed in # of grid cells radii]
        const double tau = 1.0;
        const double numResourcesToReproduce = 0.7;
        
        // predator parameters
        float pred_sr = 4.0f; // search radius [expressed in # of grid cells radii]
        int handlingTime = 5; 
        const double T = 3.0; // ammount of time that a predator can remain undetected
        const double numPreyToReproduce = 2.0; // numerical response

        // evolution parameters
        const double initLambda = 0.5;
        const double mutationRate = 0.005;
        const double stdMutation = 0.05;

        // randomness related
        unsigned int seed = 111144;

        // variables for saving and checking the results
        const int unsigned runTimeCheck = 10;
        const int unsigned saveResult = 50; 

        const int numIndividualsMonitored = 1000;

        // set outputs
        //std::string outputFileName1 = "individual_movement.csv"; 
        std::string popCyclesFileName = "population_cycles.csv";
        std::string lambdaEvolFileName = "lambda_evolution.csv";
        //std::string weightEvolFileName = "weights_evolution.csv";
        //std::string reactionNormFileName = "reaction_norm.csv";
    };


    //*** define enitites *****************************************************
    struct Prey
    {
        torus::vec_t pos = {0,0};    // position & search-radius  
        double uptake = 0.0;      // resource uptake; -1.0: dead
        double lambda = 0.0;

        Brain brain;
    };

    struct Pred
    {
        torus::vec_t pos = {0,0};    // position & search-radius  
        int catches = 0;
        int waitingTime = 0;
    };

    //*** define simulation *****************************************************
    class Simulation
    {
        using search_tree_t = torus::hrtree_t;

    public:
        Simulation(const Param&); // intiliase unknowing populations
        void run();

    private:
        void single_step(size_t);
        void random_walks();
        double graze();
        std::pair<int,int> hunt();
        void reproducePred();
        void reproducePrey();
        //void reproducePreyBrain();
        //void reproducePreyBrainFixed();
        //std::vector<int> computeVigilance();
        //void computeVigilanceRes();

        torus::grid_t<double> grid_;
        std::vector<Prey> prey_;
        std::vector<Pred> pred_;
        search_tree_t prey_tree_; 
        //search_tree_t pred_tree_;
        Param param_;

        std::mt19937_64 rng_;
        //std::ofstream ofs1_;
        std::ofstream popCycles_;
        std::ofstream lambdaEvol_;
        //std::ofstream weightEvol_;
        //std::ofstream reactionNorm_;

        int generation_;
    };
}

#endif
