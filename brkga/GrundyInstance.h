#ifndef GRUNDYINSTANCE_H
#define GRUNDYINSTANCE_H

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

class GrundyInstance {
public:
	typedef std::runtime_error Error;

	GrundyInstance(const std::string& instanceFile) throw(Error);
	virtual ~GrundyInstance();

	std::vector< std::vector<int> > nodes;
	std::vector< std::vector<int> > adjMatrix;
	std::vector<int> colorNodes;
	std::vector <std::vector<double> > initialSolutions;
	std::string instanceFileName;
	std::string solverType;
	int generation=0;
	int executionId;
	int globalBestFitness = 0;
	int firstImprovement = 0;
	int enableLoop = 0;
	int timeLimit=0;
	clock_t begin_time;
	// Getters:
	unsigned getNumNodes() const;
	void writeOutput(const std::string outputFile, const std::string instanceFile, 
		const int maxColor, const std::string solverType, bool validSolution, double executionTime,
		double pe, double pm, double rhoe, int generation, unsigned relGeneration, double relGenerationTime,
		int contRestart, int contLS);
	void readInitialSolutions(const std::string& solFile) throw(Error);

	int getMaxDegree() const;
	int getNumberOfNeighbors(int vertex) const;
	void setGlobalBestFitness(int globalBestFitness);
	void createChromossome(std::vector< int > nodesOrder, std::vector< double > &newChromossome);
	void chooseChromossome(int chooseType, int size, int p, double pe, std::vector<int> &chromosomePositions);
private:

	std::string name;
	std::string comment;
	std::string problemType;
	int maxDegree;

	unsigned nNodes;

	void readNodes(const std::string& line, int simpleInstancePattern) throw(Error);
	bool isEOF(const std::string& line) const;
	void trim(std::string& str) const;
	std::vector<double> convertToChromossome(std::vector<int> nodesOrder);
};

#endif