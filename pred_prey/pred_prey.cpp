#include <iostream>
//#include <game_watches.hpp>
#include "model.hpp"
#include <iostream>
#include "ann2.hpp"
#include <sstream>



int main()
{
	try {
		// *** Replicates for specific handling time **************************************************
		//model::Param p3;
		//for (double radious = 3.9; radious < 4.2; radious += 0.1) {
		//	p3.seed = 111141;
		//	for (int rep = 5; rep < 15; ++rep) {
		//		p3.pred_sr = radious;
		//		std::ostringstream oss;
		//		std::ostringstream oss2;
		//		std::ostringstream oss3;
		//		std::ostringstream oss4;
		//		oss << "brain" << "_rep" << rep << "_prob" << p3.mutationRate << "_mag" << p3.stdMutation << "_r" << radious << "_h" << 5 << "_population_cycles.csv";
		//		oss2 << "brain" << "_rep" << rep << "_prob" << p3.mutationRate << "_mag" << p3.stdMutation << "_r" << radious << "_h" << 5 << "_lambda_evolution.csv";
		//		oss3 << "brain" << "_rep" << rep << "_prob" << p3.mutationRate << "_mag" << p3.stdMutation << "_r" << radious << "_h" << 5 << "_weights_evolution.csv";
		//		oss4 << "brain" << "_rep" << rep << "_prob" << p3.mutationRate << "_mag" << p3.stdMutation << "_r" << radious << "_h" << 5 << "_reaction_norm.csv";
		//		p3.popCyclesFileName = oss.str();
		//		p3.lambdaEvolFileName = oss2.str(); 
		//		p3.weightEvolFileName = oss3.str();
		//		p3.reactionNormFileName = oss4.str();

		//		model::Simulation sim(p3);
		//		std::cout << "*********** Start simulation **************\n\n";
		//		sim.run();
		//		std::cout << "\n*********** End simulation ***********\n\n";
		//		oss.str(std::string()); // clean oss
		//		oss2.str(std::string()); // clean oss
		//		oss3.str(std::string()); // clean oss
		//		oss4.str(std::string()); // clean oss

		//		p3.seed += 1;
		//	}
		//}

		// *** Single simulation **********************************************************************
		model::Param p1;
		std::ostringstream oss;
		std::ostringstream oss2;
		std::ostringstream oss3;
		std::ostringstream oss4;
		oss << "evol" << "_rep" << 1 << "_prob" << p1.mutationRate << "_mag" << p1.stdMutation << "_r" << p1.pred_sr << "_h" << p1.handlingTime << "_population_cycles.csv";
		oss2 << "evol" << "_rep" << 1 << "_prob" << p1.mutationRate << "_mag" << p1.stdMutation << "_r" << p1.pred_sr << "_h" << p1.handlingTime << "_lambda_evolution.csv";
		//oss3 << "brain" << "_rep" << 5 << "_prob" << p1.mutationRate << "_mag" << p1.stdMutation << "_r" << p1.pred_sr << "_h" << 5 << "_weights_evolution.csv";
		//oss4 << "brain" << "_rep" << 5 << "_prob" << p1.mutationRate << "_mag" << p1.stdMutation << "_r" << p1.pred_sr << "_h" << 5 << "_reaction_norm.csv";
		p1.popCyclesFileName = oss.str();
		p1.lambdaEvolFileName = oss2.str(); 
		//p1.weightEvolFileName = oss3.str();
		//p1.reactionNormFileName = oss4.str();

		model::Simulation sim(p1);
		std::cout << "*********** Start simulation **************\n\n";
		sim.run();
		std::cout << "\n*********** End simulation ***********\n\n";


	}
	catch (std::exception& error) {
		std::cerr << "Error: " << error.what() << ".\n";
		exit(1);
	}
	catch (...) {
		std::cerr << "Error: unknown error occured.\n";
		exit(1);
	}
}


