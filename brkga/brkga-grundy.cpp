/**
 * The MIT License (MIT)
 *
 * Copyright (c) 2018
 * Rodrigo Franco Toso (rfrancotoso@gmail.com) and
 * Mauricio G.C. Resende
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the "Software"), to deal in
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to do
 * so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 */

#include <iostream>
#include <algorithm>
#include <string.h>
#include <vector>
#include <fstream>
#include <climits>
#include <queue>
#include <cmath>
#include "brkgaAPI/BRKGA.h"
#include "brkgaAPI/MTRand.h"

#include "GrundySolver.h"
#include "GrundyDecoder.h"
#include "GrundyInstance.h"
#include "GrundyVerifier.h"

int main(int argc, char* argv[]) {
	if(argc < 2) { std::cerr << "usage: <GrundyLIB-file>" << std::endl; return -1; }

	//std::cout << "Welcome to the BRKGA API sample driver.\nFinding a (heuristic) minimizer for "
	//		<< " the Grundy Number." << std::endl;

	std::string instanceFileName = "";
	std::string solverType = "grundyNP";
	std::string outputFile = "outputTestCon.sol";
	std::string outputFileExec = "outputTestExecCon.sol";
	double timeLimit = 300;
	double readPe = 0.30;
	double readPm = 0.1;
	double readRhoe = 0.6;
	int numThreads = 1;
	int numGenerations = INT_MAX;
	int numPopulations = 1;
	double populationFactor = 3;
	int seed = 100;
	int nChromossomeToImprove = 5;
	int randomChromossome = 2;
	int firstImprovement = 1;
	int enableReset = 0;
	int enableLoop = 0;
	int enableLS = 0;

	for(int i=0; i<argc; i++){
		if(strcmp(argv[i],"--time") == 0){
			timeLimit = atof(argv[i+1]);
			i++;
		}else if(strcmp(argv[i],"--inst") == 0){
			instanceFileName = std::string(argv[i+1]);
			i++;
		}else if(strcmp(argv[i],"--pe") == 0){
			readPe = atof(argv[i+1]);
			i++;
		}else if(strcmp(argv[i],"--pm") == 0){
			readPm = atof(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--rhoe") == 0){
			readRhoe = atof(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--solver") == 0){
			solverType = std::string(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--output") == 0){
			outputFile = std::string(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--outputExec") == 0){
			outputFileExec = std::string(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--generation") == 0){
			numGenerations = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--thread") == 0){
			numThreads = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--population") == 0){
			numPopulations = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--factor") == 0){
			populationFactor = atof(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--seed") == 0){
			seed = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--nToImprove") == 0){
			nChromossomeToImprove = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--improvementType") == 0){
			randomChromossome = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--firstImprovement") == 0){
			firstImprovement = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--reset") == 0){
			enableReset = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--loop") == 0){
			enableLoop = atoi(argv[i+1]); 
			i++;
		}else if(strcmp(argv[i],"--ls") == 0){
			enableLS = atoi(argv[i+1]); 
			i++;
		}

	}		

	const clock_t begin = clock();
	
	const std::string instanceFile = instanceFileName;
	std::cout << "Instance file: " << instanceFile << std::endl;

	if(instanceFile[instanceFile.length()-1] == 'n'){
		std::cout<<"exit\n";
		exit(1);
	}
	/*
	std::cout<<"instanceFile = "<<instanceFile<<"\n";
	std::cout<<"timeLimit = "<<timeLimit<<"\n";
	std::cout<<"pe = "<<readPe<<"\n";
	std::cout<<"pm = "<<readPm<<"\n";
	std::cout<<"rhoe = "<<readRhoe<<"\n";
	std::cout<<"solverType = "<<solverType<<"\n";
	std::cout<<"outputFile = "<<outputFile<<"\n";
	std::cout<<"generation = "<<numGenerations<<"\n";
	std::cout<<"population = "<<numPopulations<<"\n";
	std::cout<<"factor = "<<populationFactor<<"\n";
	*/

	// Read the instance:
	GrundyInstance instance(instanceFile); 	// initialize the instance
	instance.instanceFileName = instanceFile;
	instance.executionId = seed;
	instance.solverType = solverType;
	instance.begin_time = begin;
	instance.firstImprovement = firstImprovement;
	instance.enableLoop = enableLoop;
	instance.timeLimit = timeLimit;
	/*
	std::cout << "Instance read; here's the info:"
			<< "\n\tName: " << instance.getName()
			<< "\n\tComment: " << instance.getComment()
			<< "\n\tDimension: " << instance.getNumNodes()
			<< "\n\tEdge type: " << instance.getProblemType()
			<< "\n\tEdge Weight Type: " << instance.getEdgeWeightType() << std::endl;
	*/

	std::string decSolverType = solverType.compare("grundyCPopt") == 0 ? "grundyCP" : solverType;
	decSolverType = solverType.compare("grundyNPopt") == 0 ? "grundyNP" : decSolverType;
	GrundyDecoder decoder(instance, decSolverType);		// initialize the decoder

	const long unsigned rngSeed = seed;	// seed to the random number generator
	MTRand rng(rngSeed);					// initialize the random number generator

	const unsigned n = instance.getNumNodes();// size of chromosomes
	int populationSizeC = (int) populationFactor*n;
	const unsigned p = std::max(100, populationSizeC);		// size of population
	const double pe = readPe;		// fraction of population to be the elite-set
	const double pm = readPm;		// fraction of population to be replaced by mutants
	const double rhoe = readRhoe;	// probability that offspring inherit an allele from elite parent
	const unsigned K = numPopulations;		// number of independent populations
	const unsigned MAXT = numThreads;	// number of threads for parallel decoding
	
	// initialize the BRKGA-based heuristic
	
	BRKGA< GrundyDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

	// BRKGA inner loop (evolution) configuration: Exchange top individuals
	const unsigned X_INTVL = INT_MAX;	// exchange best individuals at every 100 generations
	const unsigned X_NUMBER = 2;	// exchange top 2 best
	const unsigned MAX_GENS = numGenerations;	// run for 1000 gens

	// BRKGA evolution configuration: restart strategy
	unsigned relevantGeneration = 1;	// last relevant generation: best updated or reset called
	const unsigned RESET_AFTER = 2000;
	std::vector< double > bestChromosome;
	double bestFitness = 0;
	
	std::vector<unsigned> relGenerations;
	std::vector<double> relGenerationsTime;
	int contRestart = 0;
	int contLS = 0;
	// Print info about multi-threading:
	#ifdef _OPENMP
		std::cout << "Running for " << MAX_GENS << " generations using " << MAXT
				<< " out of " << omp_get_max_threads()
				<< " available thread units..." << std::endl;
	#endif
	#ifndef _OPENMP
		std::cout << "Running for " << MAX_GENS
				<< " generations without multi-threading..." << std::endl;
	#endif
	
	
	// Run the evolution loop:
	double execTime;
	unsigned generation = 1;		// current generation
	std::ofstream ofs;
	//ofs.open(outputFileExec.c_str(), std::ofstream::out | std::ofstream::app);
	do {
		instance.generation = generation;
		algorithm.evolve();	// evolve the population for one generation
		
		// Bookeeping: has the best solution thus far improved?
		if(algorithm.getBestFitness() > bestFitness) {
			// Save the best solution to be used after the evolution chain:
			relevantGeneration = generation;
			bestFitness = algorithm.getBestFitness();
			instance.setGlobalBestFitness(bestFitness);
			bestChromosome = algorithm.getBestChromosome();
			
			if(enableLS == 1 && (solverType.compare("grundyCPopt") == 0 || solverType.compare("grundyNPopt") == 0)){
				contLS++;
				std::vector<int> chromosomePositions;
				instance.chooseChromossome(randomChromossome, nChromossomeToImprove, p, pe, chromosomePositions);
				
				for(int i=0; i<chromosomePositions.size(); i++){
					
					std::vector< double > actualChromosome = algorithm.getPopulation(0).getChromosome(chromosomePositions[i]);
					GrundySolver searchSolution(instance, actualChromosome, solverType);
					
					std::vector< double > newChromossome;
					newChromossome.resize(instance.getNumNodes());
					
					instance.createChromossome(searchSolution.getNodesOrder(), newChromossome);
					
					algorithm.getPopulation(0).setChromosome(p-1-i, newChromossome, searchSolution.getHighestColor());
				}
				algorithm.getPopulation(0).sortFitnessAfterInsert();
				bestFitness = algorithm.getBestFitness();
				instance.setGlobalBestFitness(bestFitness);
			}
			
			clock_t actual_time = clock();
			execTime = ((double) (actual_time - begin)) / (CLOCKS_PER_SEC);
			relGenerations.push_back(relevantGeneration);
			relGenerationsTime.push_back(execTime);
			/*
			std::cout << "\t" << generation
					<< ") Improved best solution thus far: "
					<< bestFitness << std::endl;
			*/
		}

		//  Evolution strategy: restart

		if(generation - relevantGeneration > RESET_AFTER && enableReset == 1) {
			//std::cout<<"reset\n";
			contRestart++;
			std::vector< double > actualChromosome = algorithm.getPopulation(0).getChromosome(0);
			GrundySolver bestOldSolution(instance, actualChromosome, decSolverType);
			algorithm.reset();	// restart the algorithm with random keys
			relevantGeneration = generation;
			algorithm.getPopulation(0).setChromosome(0, actualChromosome, bestOldSolution.getHighestColor());
			algorithm.getPopulation(0).sortFitnessAfterInsert();
			bestFitness = algorithm.getBestFitness();
			instance.setGlobalBestFitness(bestFitness);
			
			/*
			std::cout << "\t" << generation << ") Reset at generation "
					<< generation << std::endl;
			*/
		}

		// Evolution strategy: exchange top individuals among the populations
		if(generation % X_INTVL == 0 && relevantGeneration != generation) {
			//algorithm.exchangeElite(X_NUMBER);
			
			/*
			std::cout << "\t" << generation
					<< ") Exchanged top individuals." << std::endl;
			*/
		}

		// Next generation?
		//std::cout<<"generation "<<generation<<" "<<bestFitness<<std::endl;
		//ofs << instanceFile << ";" << generation << ";"	<< bestFitness << ";"	<< seed << ";" << execTime << "\n";

		++generation;
		clock_t actual_time = clock();
		execTime = ((double) (actual_time - begin)) / (CLOCKS_PER_SEC);
		//std::cout<<"execTime:"<<execTime<<"|timeLimit:"<<timeLimit<<"\n";
	} while (generation < MAX_GENS && (execTime < timeLimit));
	//ofs.close();
	// print the fitness of the top 10 individuals of each population:
	/*
	std::cout << "Fitness of the top 10 individuals of each population:" << std::endl;
	const unsigned bound = std::min(p, unsigned(10));	// makes sure we have 10 individuals
	for(unsigned i = 0; i < K; ++i) {
		std::cout << "Population #" << i << ":" << std::endl;
		for(unsigned j = 0; j < bound; ++j) {
			std::cout << "\t" << j << ") "
					<< algorithm.getPopulation(i).getFitness(j) << std::endl;
		}
	}
	*/
	const clock_t end = clock();

	// rebuild the best solution:
	bestChromosome = algorithm.getBestChromosome();
	GrundySolver bestSolution(instance, bestChromosome, decSolverType);

	// print its distance:
	std::cout << "Best solution found has objective value = "
	 		<< bestSolution.getHighestColor() << std::endl;

	std::vector<int> nodesColor = bestSolution.getNodesColor();
	std::vector<int> nodesOrder = bestSolution.getNodesOrder();

	bool validSolutionTest = false;
	if(solverType.compare("grundyNP") == 0 || solverType.compare("grundyNPopt") == 0){
		validSolutionTest = grundySolutionVerifier(instance, nodesOrder, nodesColor);
	}else {
		validSolutionTest = grundyConnectedSolutionVerifier(instance, nodesOrder, nodesColor);
	}
	
	if(validSolutionTest){
		std::cout<<"Valid solution\n";
	}else{
		std::cout<<"Invalid solution\n";
	}

	std::cout << "BRKGA run finished in " << (end - begin) / double(CLOCKS_PER_SEC) << " s." << std::endl;
	
	double executionTime = (end - begin) / double(CLOCKS_PER_SEC);

	
	instance.writeOutput(outputFile, instanceFile, bestSolution.getHighestColor(), solverType, 
		validSolutionTest, executionTime, pe, pm, rhoe, generation, 
		relGenerations[relGenerations.size()-1], relGenerationsTime[relGenerationsTime.size()-1], contRestart, contLS);

	return 0;
}