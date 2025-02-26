#ifndef GRUNDYSOLVER_H
#define GRUNDYSOLVER_H

#include <list>
#include <limits>
#include <vector>
#include <algorithm>
#include "GrundyInstance.h"


class GrundySolver {
public:
	// The constructor 'solves' the problem in O(n log n) by transforming the chromosome into
	// a tour (by getting a permutation out of the chromosome):
	GrundySolver(const GrundyInstance& instance, const std::vector< double >& chromosome, std::string solverType);
	virtual ~GrundySolver();

	int GrundySolverNormal(const GrundyInstance& instance, const std::vector< double >& chromosome, 
		std::vector<int>& nodesColor);
	int GrundySolverConnected(const GrundyInstance& instance, const std::vector< double >& chromosome,
		std::vector<int>& nodesColor);
	int GrundySolverNormalPriority(const GrundyInstance& instance, const 
		std::vector< double >& chromosome, std::vector<int>& nodesColor);
	int GrundySolverConnectedPriority(const GrundyInstance& instance, 
		const std::vector< double >& chromosome, std::vector<int>& nodesColor);

	int GrundySolverNormalPriorityBins(const GrundyInstance& instance, const 
		std::vector< double >& chromosome, std::vector<int>& nodesColor);

	int GrundySolverConnectedWithIG(const GrundyInstance& instance, 
		const std::vector< double >& chromosome, std::vector<int>& nodesColor);

	int GrundySolverWithNeighborhood(const GrundyInstance& instance, 
		const std::vector< double >& chromosome, std::vector<int>& nodesColor, std::string solverType);	
	
	int iteratedGreedy(const GrundyInstance& instance, const std::vector< double >& chromosome, 
	std::vector<int>& colorNodes, std::vector<int> &nodesOrder, int bestColor);

	int GrundySolverWithNeighborhoodCompact(const GrundyInstance& instance, 
		const std::vector< double >& chromosome, std::vector<int>& nodesColor, std::string solverType);

	void determineFirstNeighbor(const GrundyInstance& instance, const std::vector< int > &nodesOrder, 
		std::vector< std::vector < int > >& firstNeighbor);
	
	int colorSequence(const GrundyInstance& instance, std::list<int> &sequence, std::vector<int> &colors, std::vector<int> &colorNodes,
		int positionV, int newPosition, int fPosition);	

	unsigned getHighestColor() const;		// Returns the highest color
	int chooseNodeColor(const GrundyInstance& instance, const std::vector< int >& colorNodes, int vertex) const;
	std::vector<int> getNodesColor() const;
	std::vector<int> getNodesOrder() const;
	std::vector< std::vector<int> > solutions;

private:
	unsigned color;
	std::vector<int> nodesColor;
	std::vector<int> nodesOrder;
};

#endif