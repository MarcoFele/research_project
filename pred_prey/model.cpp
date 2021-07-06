#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include "model.hpp"
#include <fstream>
#include <string>
#include <iostream>
#include <random>

using namespace torus;

namespace model {

    // initialise unknowing simulation
    Simulation::Simulation(const Param& param) :
        grid_(param.GS, param.initResourcesDensity),
        prey_(param.N_prey),
        pred_(param.N_pred),
        param_(param),
        rng_(param.seed),
        generation_(0)
    {
        std::cout << "Initialising conditions ... " << "\n\n";

        auto pdist = std::uniform_real_distribution<float>(0.0, 1.0);

        for (auto& prey : prey_) {
            prey.pos = { pdist(rng_), pdist(rng_) };
            prey.lambda = param.initLambda;
            // initial lambda 
            for (auto& weight : prey.brain) {
                weight = param.initLambda;
                break;
            }
        }

        for (auto& pred : pred_) {
            pred.pos = { pdist(rng_), pdist(rng_) };
        }

        //ofs1_.open(param.outputFileName1.c_str()); // individual movement
        //if (!ofs1_.is_open())
        //    throw std::logic_error("unable to open file.");
        //ofs1_ << "time" << "," << "individual" << "," << "type" << "," << "x" << ","  << "y" << "\n";

        popCycles_.open(param.popCyclesFileName.c_str()); // population cycles
        if (!popCycles_.is_open())
           throw std::logic_error("unable to open file.");
        popCycles_ << "seed" << "," << "initPrey" << "," << "initPred" << "," << "tau" << "," << "bigT" << "," // vigilance calculation parameters
             << "handlingTime" << "," << "predSearchRadious" << ","// pred parameters
             << "preySearchRadious" << "," << "mutationRate" << "," << "stdMutation" << ","// prey parameters and evolution
             << "time" << "," << "generation" << "," << "resources" << "," 
            << "preyNumber" << "," << "predNumber" << "," << "meanLambda" << "," << "meanSat" << "," << "captures" << "," << "inPredSight" << "," 
            << "totPreyEnc" << "," << "totPredEnc" << "\n";
    
        lambdaEvol_.open(param.lambdaEvolFileName.c_str()); // lambda evolution
        if (!lambdaEvol_.is_open())
            throw std::logic_error("unable to open file.");
        lambdaEvol_ << "seed" << "," << "initPrey" << "," << "initPred" << "," << "tau" << "," << "bigT" << "," // vigilance calculation parameters
            << "handlingTime" << "," << "predSearchRadious" << ","// pred parameters
           << "preySearchRadious" << "," << "mutationRate" << "," << "stdMutation" << ","// prey parameters and evolution
           << "time" << "," << "generation" << "," << "individual" << "," 
            << "lambda" << "," << "preyEncounters" << "," << "predEncounters" << "\n";

        //weightEvol_.open(param.weightEvolFileName.c_str()); // weights evolution
        //if (!weightEvol_.is_open())
        //    throw std::logic_error("unable to open file.");
        //weightEvol_ << "seed" << "," << "initPrey" << "," << "initPred" << "," << "tau" << "," << "bigT" << "," // vigilance calculation parameters
        //    << "handlingTime" << "," << "predSearchRadious" << ","// pred parameters
        //    << "preySearchRadious" << "," << "mutationRate" << "," << "stdMutation" << ","// prey parameters and evolution
        //    << "time" << "," << "generationNumber" << "," << "individual" << "," << "weightNumber " << "," << "weight" << "\n";

        //reactionNorm_.open(param.reactionNormFileName.c_str()); // reaction norm evolution
        //if (!reactionNorm_.is_open())
        //    throw std::exception("unable to open reaction norm output");
        //reactionNorm_ << "seed" << "," << "initPrey" << "," << "initPred" << "," << "tau" << "," << "bigT" << "," // vigilance calculation parameters
        //    << "handlingTime" << "," << "predSearchRadious" << "," // pred parameters
        //    << "preySearchRadious" << "," << "mutationRate" << "," << "sdMutation" << "," // prey parameters
        //    << "generationTime" << "," << "maxTime" << "," // general simulation parameters
        //    << "individual" << "," << "numPred" << "," << "numCos" << "," << "numRes" << "," << "lambda" << "\n";
    }

    void Simulation::run()
    {
        for (size_t t = 0; t < param_.maxTime; ++t) {
            single_step(t);

            // if run multiple simulations
            if (prey_.size() == 0 || pred_.size() == 0) 
                break;
        }

        // plot reaction norm of prey subset 
        //std::cout << "\n*********** Saving the survivor reaction norm **************\n\n";
        //double maxRes = param_.resourceDensity;
        //for (int individual = 0; individual < std::min(param_.numIndividualsMonitored, static_cast<int>(prey_.size())); ++individual) {  // reaction norm
        //    if(!(individual % param_.saveResult))
        //        std::cout << "Individual:\t" << individual << "\n";
        //    for (int numPred = 0; numPred < 10; ++numPred) {
        //        for (int numCons = 0; numCons < 10; ++numCons) {
        //            for (double res = 0.0; res < maxRes; res += maxRes / 10) {
        //                reactionNorm_ << param_.seed << "," << param_.N_prey << "," << param_.N_pred << "," << param_.tau << "," << param_.T << "," // vigilance calculation parameters
        //                    << param_.handlingTime << "," << param_.pred_sr << "," // pred parameters
        //                    << param_.prey_sr << "," << param_.mutationRate << "," << param_.stdMutation << "," // prey parameters
        //                    << param_.generationTime << "," << param_.maxTime << "," // simulation parameters
        //                    << individual << "," << numPred << "," << numCons << "," << res << "," << prey_[individual].brain(numCons, numPred)[0] << "\n"; // data
        //            }
        //        }
        //    }
        //}

        //ofs1_.close();
        popCycles_.close();
        lambdaEvol_.close();
        //weightEvol_.close();
        //reactionNorm_.close();
    }


    void Simulation::single_step(size_t t)
    {        
        // move individuals
        random_walks();

        // shuffle individuals
        shuffle(pred_.begin(), pred_.end(), rng_); 
        shuffle(prey_.begin(), prey_.end(), rng_); 

        // construct tree
        //pred_tree_.build(pred_.cbegin(), pred_.cend(), [&](auto& pred) { return aabb_t{ pred.pos, {0, 0} }; });
        prey_tree_.build(prey_.cbegin(), prey_.cend(), [&](auto& prey) { return aabb_t{ prey.pos, {0, 0} }; });

        // determine vigilance level
        //std::vector<int> vecEncounters = computeVigilance(); // first prey, second predators

        // eat resoruces
        double numResourcesEaten = graze(); 

        // regenerate resources
        for (auto& c : grid_) 
            c = std::min(param_.resourceDensity, c + param_.grow_back);

        // eat prey
        std::pair<int,int> numPreyEaten = hunt();

        // reproduce individual
        if (!(t % param_.generationTime) && t != 0) {
           reproducePrey(); 
           reproducePred();
            ++generation_;
        }

        // runtimecheck
        if (!(t % param_.runTimeCheck)) {
            double numResources = 0;
            for (auto& c : grid_) // count resources
                numResources += c;
            std::cout << "time:\t" << t
                << "\ngeneration:\t" << generation_
                << "\nresources:\t" << numResources
                << "\nprey:\t" << prey_.size()
                << "\npredators:\t" << pred_.size() 
                << "\nper capita # of prey eaten:\t" << static_cast<double>(numPreyEaten.second) / static_cast<double>(pred_.size())
                << "\napproximate per capita # resources eaten:\t" << numResourcesEaten / static_cast<double>(prey_.size()) << "\n\n";
        }
        // write output
        if (!(t % param_.saveResult)) {
            //int preyInd = 0; // individual movement
            //for_each(prey_.cbegin(), prey_.cend(), [&](auto& individual) mutable {
            //    ofs1_ << param_.handlingTime << "," << param_.pred_sr << "," // pred parameters
            //       << param_.prey_sr << "," << param_.mutationRate << "," << param_.stdMutation << "," // prey parameterst 
            //       << t << "," << preyInd++ << "," << "prey" << "," << individual.pos_sr.center << "\n";
            //    });
            //int predInd = 0;
            //for_each(pred_.cbegin(), pred_.cend(), [&](auto& individual) mutable {
            //    ofs1_ <<  << param_.handlingTime << "," << param_.pred_sr << "," // pred parameters
            //      << param_.prey_sr << "," << param_.mutationRate << "," << param_.stdMutation << "," // prey parameters
            //      << t << "," << predInd++ << "," << "pred" << "," << individual.pos << "\n";
            //    });

            double numResources = 0; // population cycles
            for (auto& c : grid_) // calculate res
                numResources += c;
            double totLambda = 0; // calculate lambda
            for (auto& prey : prey_)
                totLambda += prey.lambda;
            double totSaturation = 0; // calculate saturation
            for (auto& pred : pred_)
                totSaturation += pred.waitingTime;
            popCycles_ << param_.seed << "," << param_.N_prey << "," << param_.N_pred << "," << param_.tau << "," << param_.T << "," // vigilance calculation parameters
                << param_.handlingTime << "," << param_.pred_sr << "," // pred parameters
                << param_.prey_sr << "," << param_.mutationRate << "," << param_.stdMutation << "," // prey parameters
                << t << "," << generation_ << "," << numResources << "," << prey_.size() << "," << pred_.size() << ","
                << totLambda / static_cast<double>(prey_.size()) << ","
                << totSaturation / static_cast<double>(pred_.size()) << ","
                << numPreyEaten.second << "," // prey captures
                << numPreyEaten.first << "\n"; // << "," // in predator sights
                //<< vecEncounters[vecEncounters.size() - 2] << "," // tot prey enocunters from prey perspective
                //<< vecEncounters[vecEncounters.size() - 1] << "\n"; // tot pred enocunters from predator perspective

            int numMonitored = std::min(param_.numIndividualsMonitored, static_cast<int>(prey_.size())); // lambda evolution
            for (int i = 0; i < numMonitored; ++i) { // for prey
                lambdaEvol_ << param_.seed << "," << param_.N_prey << "," << param_.N_pred << "," << param_.tau << "," << param_.T << "," // vigilance calculation parameters
                    << param_.handlingTime << "," << param_.pred_sr << "," << param_.prey_sr << "," // pred parameters
                    << param_.mutationRate << "," << param_.stdMutation << "," // prey parameters
                    << t << "," << generation_ << "," << i << ","
                    << prey_[i].lambda << "\n"; // << ","           
                    //<< vecEncounters[i] << "," // save prey encounters
                    //<< vecEncounters[i + numMonitored] << "\n"; // save pred encounters
            }

            //for (int i = 0; i < std::min(param_.numIndividualsMonitored, static_cast<int>(prey_.size())); ++i) { // weight evolution
            //    int weightNum = 0;
            //    for (auto& weight : prey_[i].brain) {
            //        weightEvol_ << param_.seed << "," << param_.N_prey << "," << param_.N_pred << "," << param_.tau << "," << param_.T << "," // vigilance calculation parameters
            //            << param_.handlingTime << "," << param_.pred_sr << "," // pred parameters
            //            << param_.prey_sr << "," << param_.mutationRate << "," << param_.stdMutation << "," // prey parameters
            //            << t << "," << generation_ << "," << i << "," << weightNum << "," << weight << "\n";
            //        ++weightNum;
            //    }
            //}
        }
    }

    void Simulation::random_walks()
    {
        std::exponential_distribution<float> movmentLength(param_.randomWalk);
        auto adist = std::uniform_real_distribution<>(0.0, 2.0 * M_PI);
        auto sdist = [&]() {
            const auto a = adist(rng_);
            return vec_t{ static_cast<float>(std::cos(a)), static_cast<float>(std::sin(a)) };
        };
        for (auto& prey : prey_) {
            float step = movmentLength(rng_);
            while (isinf(step)) 
                step = movmentLength(rng_);
            prey.pos = wrap(prey.pos + step * 2.0f * grid_.pixel_radii() * sdist()); 
        }
        for (auto& pred : pred_) {
            float step = movmentLength(rng_);
            while(isinf(step)) 
                step = movmentLength(rng_);
            pred.pos = wrap(pred.pos + step * 2.0f * grid_.pixel_radii() * sdist());
        }
    }


    // ******** functions for prey ****************************************
    //std::vector<int> Simulation::computeVigilance() {
    //    int numMonitored = std::min(param_.numIndividualsMonitored, static_cast<int>(prey_.size())); // lambda evolution
    //    std::vector<int> enocuntersCount((numMonitored * 2) + 2, 0); // index:0-999 are individual monitored prey enc, index:1000-1999 are individual monitored pred enc, index:2000 is total ecounters prey, index:2001 is total encounters predator, in total 2002 slots (if numIndividualsMonitored == 1000)
    //    int individual = 0;
    //    auto sr = param_.prey_sr * grid_.pixel_radius();
    //    for (auto& prey : prey_) {
    //        // manage saving individual 
    //        bool saveIndividual = 0;
    //        if (individual < numMonitored)
    //            saveIndividual = 1;

    //        auto pos = prey.pos;
    //        // count neighbouring predators
    //        double numPred = 0.0;
    //        pred_tree_.query({ pos, {sr,sr} }, [&](size_t j) {
    //            auto dd = distance2(pos, pred_[j].pos);
    //            if (dd < sr * sr) {
    //                ++numPred;
    //                // count pred enocunters
    //                if (saveIndividual) {
    //                    ++enocuntersCount[individual + numMonitored]; // count predator encounters for specific individual
    //                    ++enocuntersCount[enocuntersCount.size() - 1]; // count total predator encounters
    //                }
    //            }
    //            });
    //        // count neighbouring cospecifics
    //        double numCons = 0.0;
    //        prey_tree_.query({ pos, {sr,sr} }, [&](size_t j) {
    //            auto dd = distance2(pos, prey_[j].pos);
    //            if (dd < sr * sr) {
    //                ++numCons;
    //                // count prey enocunters
    //                if (saveIndividual) {
    //                    ++enocuntersCount[individual]; // count prey encounters for specific individual
    //                    ++enocuntersCount[enocuntersCount.size() - 2]; // count total prey encounters
    //                }
    //            }
    //            });
    //        // update lambda value
    //        prey.lambda = prey.brain(std::array<double, 2>{numCons, numPred})[0];
    //        if (prey.lambda < 0)
    //            throw std::exception("negative lambda value");

    //        if (saveIndividual) 
    //            ++individual;
    //    }

    //    return enocuntersCount;
    //}

    //void Simulation::computeVigilanceRes() {
    //    auto sr = param_.prey_sr * grid_.pixel_radius();
    //    for (auto& prey : prey_) {
    //        auto pos = prey.pos;
    //        // count neighbouring predators
    //        double numPred = 0.0;
    //        pred_tree_.query({ pos, {sr,sr} }, [&](size_t j) {
    //            auto dd = distance2(pos, pred_[j].pos);
    //            if (dd < sr * sr)
    //                ++numPred;
    //            });
    //        // count neighbouring cospecifics
    //        double numCons = 0;
    //        prey_tree_.query({ pos, {sr,sr} }, [&](size_t j) {
    //            auto dd = distance2(pos, prey_[j].pos);
    //            if (dd < sr * sr)
    //                ++numCons;
    //            });
    //        // update lambda value
    //        prey.lambda = prey.brain(std::array<double, 3>{numCons, numPred, grid_(pos)})[0];
    //        if (prey.lambda < 0)
    //            throw std::exception("negative lambda value");
    //    }
    //}

    double Simulation::graze()
    {
        double numResourcesEaten = 0.0;

        // NO fair share between prey on same cell, no frequency depended selection
        for (auto& prey : prey_) {
            double feedingRate = 1.0 / (1.0 + prey.lambda * param_.tau); 
            double resourcesFound = feedingRate * grid_(prey.pos);
            grid_(prey.pos) -= resourcesFound;
            prey.uptake += resourcesFound;
            numResourcesEaten += resourcesFound;

            //## needs to be enforced
            // exception check
            if (grid_(prey.pos) < 0) {
                std::ostringstream oss;
                oss << "resources went negative in pixel x=" << prey.pos[0] << ", y=" << prey.pos[1];
                throw std::runtime_error(oss.str().c_str());
            }
        }
        return numResourcesEaten;
    }

    void Simulation::reproducePrey()
    {
        auto pdist = std::uniform_real_distribution<float> (0.0, 1.0);
        std::bernoulli_distribution mutate(param_.mutationRate);
        std::vector<Prey> newGen;

        for (const Prey& prey : prey_) {
            int count = 0;
            if (prey.uptake > 0) { // poisson distribution does not accept zero as mean
                std::poisson_distribution<> howMany(prey.uptake / param_.numResourcesToReproduce);
                count = howMany(rng_);
            }
            while (count) {
                double newLambda = prey.lambda;
                // mutation occurs
                if (mutate(rng_)) {
                    std::normal_distribution<> magnitude(prey.lambda, param_.stdMutation);
                    newLambda = std::max(magnitude(rng_), 0.0); // lambda cannot have negative values 
                }
                Prey newPrey;
                newPrey.pos = { pdist(rng_), pdist(rng_) };
                newPrey.lambda = newLambda;
                newGen.push_back(newPrey);
                --count;
            }
        }
        // runtime check and create offsprings
        //if (newGen.empty())
          // throw std::runtime_error("prey population has no offsprings");
        prey_.swap(newGen); 
    }

    //void Simulation::reproducePreyBrain() // not keep constant population size
    //{
    //    auto pdist = std::uniform_real_distribution<float>(0.0, 1.0);
    //    std::bernoulli_distribution mutate(param_.mutationRate);
    //    std::vector<Prey> newGen;

    //    for (const Prey& prey : prey_) {
    //        double preyLambda = prey.lambda;
    //        std::normal_distribution<> magnitude(0.0, (prey.lambda == 0 ? param_.stdMutation : param_.stdMutation * preyLambda));

    //        int count = 0;
    //        if (prey.uptake > 0) { // poisson distribution does not accept zero as mean
    //            std::poisson_distribution<> howMany(prey.uptake / param_.numResourcesToReproduce);
    //            count = howMany(rng_);
    //        }
    //        while (count) {
    //            Prey newPrey = prey;
    //            newPrey.uptake = 0;
    //            for (auto& weight : newPrey.brain) {
    //                if (mutate(rng_)) {
    //                    weight += magnitude(rng_);
    //                }
    //            }
    //            newPrey.pos = { pdist(rng_), pdist(rng_) };
    //            newGen.push_back(newPrey);
    //            --count;
    //        }
    //    }
    //    // runtime check and create offsprings
    //    //if (newGen.empty())
    //        //throw std::runtime_error("prey population has no offsprings");
    //    prey_.swap(newGen);
    //}

    //void Simulation::reproducePreyBrainFixed() // keep constant population size
    //{
    //    auto pdist = std::uniform_real_distribution<float>(0.0, 1.0);
    //    std::bernoulli_distribution mutate(param_.mutationRate);
    //    std::vector<Prey> newGen;
    //    std::normal_distribution<> magnitude(0.0, param_.stdMutation);
    //    std::vector<double> reproduceChance(prey_.size());
    //    for (int i = 0; i < prey_.size(); ++i) {
    //        reproduceChance[i] = prey_[i].uptake;
    //    }
    //    std::discrete_distribution<int> weightedLottery(reproduceChance.begin(), reproduceChance.end());

    //    for (int i = 0; i < param_.N_prey; ++i) {
    //        Prey newPrey = prey_[weightedLottery(rng_)];
    //        newPrey.pos = { pdist(rng_), pdist(rng_) };
    //        newPrey.uptake = 0.0;
    //        for (auto& weight : newPrey.brain) {
    //            if (mutate(rng_)) {
    //                weight += magnitude(rng_);
    //            }
    //        }
    //        newGen.push_back(newPrey);
    //    }
    //    // runtime check and create offsprings
    //    if (newGen.empty())
    //        throw std::runtime_error("prey population has no offsprings");
    //    else
    //        prey_.swap(newGen);
    //}


    // ******** functions for predator ****************************************

    std::pair<int,int> Simulation::hunt() // first encounters, second captures
    {
        int numPreyEaten = 0;
        int numPreyEnc = 0;
        auto scaledRadious = grid_.pixel_radius() * param_.pred_sr;

        std::vector<size_t> candidates;
        for (auto& pred : pred_) {
            if (!pred.waitingTime) {
                const auto pos = pred.pos;
                candidates.clear();
                prey_tree_.query({ pos, {scaledRadious, scaledRadious} }, [&](size_t j) {
                    const auto dd = distance2(pos, prey_[j].pos);
                    if (dd < scaledRadious * scaledRadious) { // if prey is close
                        if (prey_[j].uptake >= 0) {// if prey is alive
                            candidates.push_back(j);
                            ++numPreyEnc;
                        }
                    }
                    });
                shuffle(candidates.begin(), candidates.end(), rng_); // chose randomly the prey to attack
                for (std::vector<size_t>::iterator it = candidates.begin(); it != candidates.end(); ++it) {
                    double feedingRate = 1.0 / (1.0 + prey_[*it].lambda * param_.tau);
                    double deadProb = feedingRate * exp(-(param_.T * prey_[*it].lambda));
                    std::bernoulli_distribution attack(deadProb);
                    if (attack(rng_)) {
                        ++pred.catches;
                        pred.waitingTime = param_.handlingTime;
                        prey_[*it].uptake = -1.0;   // doomed
                        ++numPreyEaten;
                        break;
                    }
                }
            }
            else {
                --pred.waitingTime;
                // runtime check
                if (pred.waitingTime < 0)
                    throw std::runtime_error("impossible negative waiting time");
            }
        }
        // remove dead prey
        prey_.erase( // remove_if preserves the order but slower
            std::partition(prey_.begin(), prey_.end(), [](const auto& prey) { return prey.uptake >= 0.0; }),
            prey_.end()
        );

        //if (prey_.empty())
         // throw std::exception("prey population got totally eaten");

        return std::pair<int, int> {numPreyEnc, numPreyEaten};
    }
            
    void Simulation::reproducePred()
    {
        auto pdist = std::uniform_real_distribution<float>(0.0, 1.0);
        std::vector<Pred> newGen;
        for (Pred& pred : pred_) {
            int count = 0;
            if (pred.catches > 0) {
                std::poisson_distribution<int> howMany(pred.catches / param_.numPreyToReproduce);
                count = howMany(rng_);
            }
            while (count) {
                Pred newPred;
                newPred.pos = { pdist(rng_), pdist(rng_) };
                newGen.push_back(newPred); // generate new individual
                --count;
            }
        }
        // runtime check
        //if (newGen.empty()) 
            //throw std::runtime_error("predator population unable to reproduce");

        pred_.swap(newGen);
    }
}
