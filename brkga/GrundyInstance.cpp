#include <string.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <cmath>
#include "GrundyInstance.h"

GrundyInstance::GrundyInstance(const std::string& instanceFile) throw(GrundyInstance::Error) :
		name(""), comment(""), problemType(""), nNodes(0) {
	std::ifstream fin(instanceFile.c_str());
	if(! fin) { throw Error("GrundyInstance: Cannot open input file."); }

	std::string line;
	std::string token;
	int simpleInstancePattern = 0;
	int nEdges=0;
	int auxNumNodes, auxNumEdges;
	// WARNING: the code below assumes an ordered input file following a280.tsp
	//          it will not work with all instances in the TSPLIB.
	try {
		if(strstr(instanceFile.c_str(), "B_") != NULL || strstr(instanceFile.c_str(), "CB_") != NULL){
			std::getline(fin, line);
			std::istringstream sin(line);
			sin >> nNodes; 
			simpleInstancePattern = 1;
		}else if(strstr(instanceFile.c_str(), "mtx") != NULL){
			std::getline(fin, line);
			std::getline(fin, line);
			std::istringstream sin(line);
			sin >> nNodes >> auxNumNodes >> auxNumEdges; 
			simpleInstancePattern = 2;
		}else{
			char beginCharRead = 'p', actualCharRead = 'c';
			std::string temp = "";

			while(actualCharRead != beginCharRead){
				std::getline(fin, line);
				actualCharRead = line[0];
			}
			std::istringstream sin(line);
			sin >> beginCharRead >> temp >> nNodes >> nEdges;
		}
		
		nodes.resize(nNodes);
		colorNodes.resize(nNodes);
		adjMatrix.resize(nNodes);
		
		std::getline(fin, line);
		while(!isEOF(line) && !fin.eof()) {
			readNodes(line, simpleInstancePattern);
			std::getline(fin, line);
		}

		for(int i=0; i<nNodes; i++){
			adjMatrix[i].resize(nNodes,0);
			adjMatrix[i][i] = 1;
			for(int j=0; j < nodes[i].size(); j++){
				adjMatrix[i][nodes[i][j]] = 1;
			}
		}

		maxDegree = 0;
		for(int i = 0; i < nNodes; i++){
			int vSize =  nodes[i].size();
			maxDegree = std::max(maxDegree, vSize);
		}

	}
	catch(const Error& error) { throw error; }
}

void GrundyInstance::readInitialSolutions(const std::string& solFile) 
	throw(GrundyInstance::Error) {

	std::ifstream fin(solFile.c_str());
	if(! fin) { throw Error("GrundyInstance: Cannot open input file."); }

	char identifyType;
	int readNode;
	std::string line;
	std::string token;
	try {
		std::vector<int> nodesInitialSol;
		while(!fin.eof()){
			std::getline(fin, line);
			std::istringstream sin(line);
			sin >> identifyType;
			if(identifyType == 'v'){
				sin >> readNode;
				nodesInitialSol.push_back(readNode);
			}else{
				initialSolutions.push_back(convertToChromossome(nodesInitialSol));
				nodesInitialSol.clear();
			}
		}
	}catch(const Error& error) { throw error; }
}

std::vector<double> GrundyInstance::convertToChromossome(std::vector<int> nodesOrder){
	std::vector<double> encodedChromossome;
	encodedChromossome.resize(nNodes);
	double step = 1/nNodes, alleloValue = 1.0;
	for(int i=0; i<nodesOrder.size(); i++){
		encodedChromossome[nodesOrder[i]] = alleloValue;
		alleloValue -= step;
	}
	return encodedChromossome;
}

GrundyInstance::~GrundyInstance() { }

unsigned GrundyInstance::getNumNodes() const { return nNodes; }

void GrundyInstance::setGlobalBestFitness(int newGlobalBestFitness) { globalBestFitness = newGlobalBestFitness;}

void GrundyInstance::createChromossome(std::vector< int > nodesOrder, std::vector< double > &newChromossome){
	std::vector< double > randomNumbers;
	for(int i=0; i<nNodes; i++){
		double aux = (double) (rand()%10000) / 10000;
		randomNumbers.push_back(aux);
	}
	std::sort(randomNumbers.begin(), randomNumbers.end());
	std::reverse(randomNumbers.begin(), randomNumbers.end());
	for(int i=0; i<nodesOrder.size(); i++){
		newChromossome[nodesOrder[i]] = randomNumbers[i];
	}
}

void GrundyInstance::readNodes(const std::string& line, int simpleInstancePattern) throw(GrundyInstance::Error) {
	std::istringstream sin(line);

	int x, y;
	std::string edgeId;

	if(simpleInstancePattern > 0){
		sin >> x >> y;
		if(simpleInstancePattern == 1){
			x++;
			y++;
		}
	}else{
		sin >> edgeId >> x >> y;
	}
	nodes[x-1].push_back(y-1);
	nodes[y-1].push_back(x-1);
}

bool GrundyInstance::isEOF(const std::string& line) const {
	if(line.find("EOF") == std::string::npos || line == "") { return false; }
	return true;
}

void GrundyInstance::trim(std::string& str) const {
	// trim white spaces at the beginning:
	unsigned begin = 0;
	while(begin < str.size() && str[begin] == ' ') { ++begin; }

	unsigned end = str.size() - 1;
	while(end > 0 && str[end] == ' ') { --end; }

	str = str.substr(begin, end - begin + 1);
}

void GrundyInstance::writeOutput(const std::string outputFile, const std::string instanceFile, 
	const int maxColor, const std::string solverType, bool validSolution, double executionTime, 
	double pe, double pm, double rhoe, int generation, unsigned relGeneration, double relGenerationTime,
	int contRestart, int contLS){

	std::ofstream ofs;
	ofs.open(outputFile.c_str(), std::ofstream::out | std::ofstream::app);

	ofs << instanceFile << ";" << solverType << ";" << maxColor << ";" << validSolution << ";" << executionTime << ";" 
	<< pe << ";" << pm << ";" << rhoe << ";" << generation << ";" << relGeneration << ";" 
	<< relGenerationTime << ";" << contRestart << ";" << contLS << "\n";
	ofs.close();
}

int GrundyInstance::getNumberOfNeighbors(int vertex) const {
	return nodes[vertex].size();
}

int GrundyInstance::getMaxDegree() const {
	return maxDegree;
}

void GrundyInstance::chooseChromossome(int chooseType, int size, int p, double pe, std::vector< int > &chromosomePositions) {
	int startNonElite = ((int) p * pe) + 1;
	if(chooseType == 0){
		for(int i=0; i<size; i++){
			chromosomePositions.push_back(i);
		}
	}else if(chooseType == 1){
		for(int i=0; i<ceil(size/2); i++){
			chromosomePositions.push_back(i);
		}
		for(int i=0; i<floor(size/2); i++){
			chromosomePositions.push_back(startNonElite + i);
		}
	}else if(chooseType == 2){
		std::set<int> setP;
		setP.insert(0);
		while(setP.size() < size){
			setP.insert(rand() % startNonElite);
		}
		for(std::set<int>::iterator it=setP.begin(); it!=setP.end(); it++){
			chromosomePositions.push_back(*it);
		}
	}
}