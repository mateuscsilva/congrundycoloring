#include "GrundySolver.h"
#include <iostream>
#include <cmath>
#include <set>
#include <queue>  
#include <utility>
#include <string>
#include <algorithm>
#include <climits>
#include <cmath>
#include <list>

GrundySolver::GrundySolver(const GrundyInstance& instance, const std::vector< double >& chromosome,
	std::string solverType) : color(0), nodesColor(instance.getNumNodes()) {
	
	if(solverType.compare("grundyNP") == 0){
		color = GrundySolverNormalPriority(instance, chromosome, nodesColor);
	}else if(solverType.compare("grundyCP") == 0){
		color = GrundySolverConnectedPriority(instance, chromosome, nodesColor);
	}else if(solverType.compare("grundyCPopt") == 0 || solverType.compare("grundyNPopt") == 0){
		color = GrundySolverWithNeighborhoodCompact(instance, chromosome, nodesColor, solverType);
	}
}

GrundySolver::~GrundySolver() {
}


int GrundySolver::GrundySolverNormalPriority(const GrundyInstance& instance, 
	const std::vector< double >& chromosome, std::vector<int>& nodesColor){
	
	int color = 0;
	std::pair <double, int> pairPriorityVertex;
	std::vector <int> colorNodes;
	std::priority_queue<std::pair<double, int> > notSelectedNodes;
	
	colorNodes.resize(instance.getNumNodes());
	

	for(int i = 0; i < instance.getNumNodes(); i++){ 
		std::pair <double, int> pairChromPriority;
		pairChromPriority = std::make_pair(chromosome[i], i);
		notSelectedNodes.push(pairChromPriority); 
	}
	
	while(!notSelectedNodes.empty()) {
		pairPriorityVertex = notSelectedNodes.top();
		int vertex = pairPriorityVertex.second;
		colorNodes[vertex] = chooseNodeColor(instance, colorNodes, vertex);
		notSelectedNodes.pop();	
		nodesOrder.push_back(vertex);

		if(colorNodes[vertex] > color){
			color = colorNodes[vertex];
		}
	}
	nodesColor = colorNodes;
	return color;
}


int GrundySolver::GrundySolverConnectedPriority(const GrundyInstance& instance, 
	const std::vector< double >& chromosome, std::vector<int>& nodesColor){
	
	int bestColor = 0, countColoredVertex = 0;
	std::pair <double, int> pairPriorityVertex;
	std::vector <int> colorNodes;
	std::priority_queue<std::pair<double, int> > notSelectedNodes;
	std::vector<int> alreadyEnqueued;
	
	colorNodes.resize(instance.getNumNodes());
	alreadyEnqueued.resize(instance.getNumNodes());

	while(countColoredVertex < instance.getNumNodes()){
		int color = 0;
		int maxValue = -1;
		int maxValuePos = 0;
		for(int i = 0; i < chromosome.size(); i++){ 
			if(maxValue < chromosome[i] && alreadyEnqueued[i] == 0){
				maxValuePos = i;
				maxValue = chromosome[i];	
			}
		}

		int firstVertex = maxValuePos; 
		colorNodes[firstVertex] = 1;
		alreadyEnqueued[firstVertex] = 1;
		countColoredVertex++;
		nodesOrder.push_back(firstVertex);
		
		for(int i = 0; i < instance.nodes[firstVertex].size(); i++){
			std::pair <double, int> pairChromPriority;
			pairChromPriority = std::make_pair(chromosome[instance.nodes[firstVertex][i]], 
			instance.nodes[firstVertex][i]);
			notSelectedNodes.push(pairChromPriority); 
			alreadyEnqueued[instance.nodes[firstVertex][i]] = 1;
		}

		while(!notSelectedNodes.empty()) {
			pairPriorityVertex = notSelectedNodes.top();
			int vertex = pairPriorityVertex.second;
			colorNodes[vertex] = chooseNodeColor(instance, colorNodes, vertex);
			countColoredVertex++;
			notSelectedNodes.pop();	
			nodesOrder.push_back(vertex);

			for(int j = 0; j < instance.nodes[vertex].size(); j++){
				if(alreadyEnqueued[instance.nodes[vertex][j]] == 0){
					std::pair <double, int> pairChromPriority;
					pairChromPriority = std::make_pair(chromosome[instance.nodes[vertex][j]], 
					instance.nodes[vertex][j]);
					notSelectedNodes.push(pairChromPriority); 
					alreadyEnqueued[instance.nodes[vertex][j]] = 1;
				}
			}

			if(colorNodes[vertex] > color){
				color = colorNodes[vertex];
			}
		}

		if(bestColor < color){
			bestColor = color;
		}	
	}

	nodesColor = colorNodes;
	return bestColor;
}

int GrundySolver::GrundySolverWithNeighborhoodCompact(const GrundyInstance& instance, 
		const std::vector< double >& chromosome, std::vector<int>& nodesColor, std::string solverType){
		
	std::ofstream ofs;
	int initialSol = 0;
	if(solverType.compare("grundyCPopt") == 0){
		initialSol = GrundySolverConnectedPriority(instance, chromosome, nodesColor);
	}else if(solverType.compare("grundyNPopt") == 0){
		initialSol = GrundySolverNormalPriority(instance, chromosome, nodesColor);
	}
	int saveInitalSol = initialSol;
	bool hasIncrease = true;
	int iter=0;
	std::vector<int> nodesOrderCopy;
	std::vector< int > nodesPosition;
	std::vector< int > colors;
	nodesPosition.resize(instance.getNumNodes());
	clock_t actual_time = clock();
	double	execTime = ((double) (actual_time - instance.begin_time)) / (CLOCKS_PER_SEC);
	
	nodesOrderCopy = nodesOrder;
	while(hasIncrease && (execTime < instance.timeLimit)){
		nodesOrderCopy = nodesOrder;
		std::vector< std::vector<int> > betterOrders;
		std::vector< std::vector < int > > firstNeighbor;
		int bestNewSol = -1;
		int bestPosition = -1;
		hasIncrease=false;

		if(solverType.compare("grundyCPopt") == 0){
			determineFirstNeighbor(instance, nodesOrderCopy, firstNeighbor);
		}
		
		for(int i=0; i<nodesOrderCopy.size(); i++){
			nodesPosition[nodesOrderCopy[i]]=i;
		}
		
		bool firstImprovementStop = false;
		for(int i=0; i<instance.getNumNodes(); i++){
			for(int j=0; j<instance.nodes[i].size(); j++){
				//move left
				if(nodesPosition[i] > nodesPosition[instance.nodes[i][j]]){
					
					if(solverType.compare("grundyCPopt") == 0 && 
					nodesPosition[instance.nodes[i][j]] < firstNeighbor[nodesPosition[i]][0] && 
					nodesPosition[instance.nodes[i][j]] != 0){
						continue;
					}

					int position=-1;
					nodesOrder.clear();
					std::list<int> sequence(nodesOrderCopy.begin(), nodesOrderCopy.end());
					std::vector<int> colorNodes;
					int neighborhoodSol = colorSequence(instance, sequence, nodesColor, colorNodes, nodesPosition[i], 
						nodesPosition[instance.nodes[i][j]], nodesPosition[instance.nodes[i][j]]-1);
					copy(sequence.begin(),sequence.end(),back_inserter(nodesOrder));

					if(neighborhoodSol > initialSol){
						betterOrders.push_back(nodesOrder);
						if(neighborhoodSol > bestNewSol){
							bestNewSol = neighborhoodSol;
							bestPosition = betterOrders.size()-1;
							colors=colorNodes;
						}
						if(instance.firstImprovement == 1){
							firstImprovementStop = true;
							break;
						}
					}
				}

				//move right
				if(nodesPosition[i] < nodesPosition[instance.nodes[i][j]]){

					if(solverType.compare("grundyCPopt") == 0){
						bool jumpNoConnected = false;
						for(int l=0; l<instance.nodes[i].size(); l++){
							if(nodesPosition[instance.nodes[i][l]] <= nodesPosition[instance.nodes[i][j]] &&
							nodesPosition[instance.nodes[i][l]] > nodesPosition[i] &&
							firstNeighbor[instance.nodes[i][l]][0] == nodesPosition[i] && 
							firstNeighbor[instance.nodes[i][l]].size() == 1){
								jumpNoConnected = true;
								break;
							}
						}
						if(jumpNoConnected){
							continue;
						}
					}

					int position=-1;
					nodesOrder.clear();
					std::list<int> sequence(nodesOrderCopy.begin(), nodesOrderCopy.end());
					std::vector<int> colorNodes;
					int neighborhoodSol = colorSequence(instance, sequence, nodesColor, colorNodes, nodesPosition[i], 
						nodesPosition[instance.nodes[i][j]], nodesPosition[i]-1);
					copy(sequence.begin(),sequence.end(),back_inserter(nodesOrder));
					
					if(neighborhoodSol > initialSol){
						betterOrders.push_back(nodesOrder);
						if(neighborhoodSol > bestNewSol){
							bestNewSol = neighborhoodSol;
							bestPosition = betterOrders.size()-1;
							colors=colorNodes;
						}
						if(instance.firstImprovement == 1){
							firstImprovementStop = true;
							break;
						}
					}
				}
			}
			if(firstImprovementStop){
				break;
			}
		}
		if(bestPosition == -1 || bestNewSol <= initialSol){
			hasIncrease=false;
		}else{
			initialSol = bestNewSol;
			if(instance.enableLoop == 1){
				hasIncrease=true;
			}else{
				hasIncrease=false;
			}
			nodesOrder = betterOrders[bestPosition];
			nodesColor=colors;
			nodesOrderCopy = nodesOrder;
		}
		actual_time = clock();
		execTime = ((double) (actual_time - instance.begin_time)) / (CLOCKS_PER_SEC);
		iter++;
	}
	nodesOrder = nodesOrderCopy;
	actual_time = clock();
	execTime = ((double) (actual_time - instance.begin_time)) / (CLOCKS_PER_SEC);
	return initialSol;
}

int GrundySolver::colorSequence(const GrundyInstance& instance, std::list<int> &sequence, std::vector<int> &colors, std::vector<int> &colorNodes,
	int positionV, int newPosition, int fPosition){
		colorNodes.resize(instance.getNumNodes());
		//update to new sequence
		std::list<int>::iterator it = sequence.begin();
		advance(it, positionV);
		int vertex = (*it);
		sequence.erase(it);
		it = sequence.begin();
		advance(it, newPosition);
		sequence.insert(it, vertex);

		//copy elements before critical interval
		it = sequence.begin();

		for(int i=0; i<=fPosition; i++){
			colorNodes[(*it)] = colors[(*it)];
			it++;
		}

		//check change in critical interval
		bool checkChange = false;
		int intervalBegin = 0;
		int intervalEnd = instance.getNumNodes()-1;
		if(newPosition < positionV){
			intervalBegin = newPosition;
			intervalEnd = newPosition;
		}else{
			intervalBegin = positionV;
			intervalEnd = newPosition;
		}

		for(int i=intervalBegin; i<=intervalEnd; i++){
			if(instance.adjMatrix[vertex][(*it)] == 1 || checkChange){
				colorNodes[(*it)] = chooseNodeColor(instance, colorNodes, (*it));
				if(colorNodes[(*it)] != colors[(*it)]){
					checkChange = true;
				}
			}else{
				colorNodes[(*it)] = colors[(*it)];	
			}
			it++;
		}

		//after critical interval
		if(!checkChange && intervalEnd < instance.getNumNodes()-1){
			for(int i=intervalEnd+1; i<instance.getNumNodes(); i++){
				colorNodes[(*it)] = colors[(*it)];
				it++;
			}
		}else if(checkChange && intervalEnd < instance.getNumNodes()-1){
			for(int i=intervalEnd+1; i<instance.getNumNodes(); i++){
				colorNodes[(*it)] = chooseNodeColor(instance, colorNodes, (*it));
				it++;
			}
		}

		//recovery maxcolor
		int max=0;
		for(int i=0; i<colorNodes.size(); i++){
			if(max < colorNodes[i]){
				max = colorNodes[i];
			}
		}

		return max;
}

unsigned GrundySolver::getHighestColor() const { return color; }


int GrundySolver::chooseNodeColor(const GrundyInstance& instance, const std::vector< int >& colorNodes,
 	int vertex) const {
	
	int actualColor = 0;
	bool isChoosedColor = false;
	while(!isChoosedColor){
		actualColor++;
		isChoosedColor = true;
		for(int i=0; i < instance.nodes[vertex].size(); i++){
			if(actualColor == colorNodes[instance.nodes[vertex][i]]){
				isChoosedColor = false;
				break;
			}
		}
	}
	return actualColor; 
}

void GrundySolver::determineFirstNeighbor(const GrundyInstance& instance, const std::vector< int > &nodesOrder, 
	std::vector< std::vector < int > >& firstNeighbor) {
	
	firstNeighbor.resize(nodesOrder.size());
	for(int i=0; i<nodesOrder.size(); i++){
		if(firstNeighbor[nodesOrder[i]].size() == 0){
			firstNeighbor[nodesOrder[i]].push_back(i);
		}
		for(int j=0; j<instance.nodes[nodesOrder[i]].size(); j++){
			firstNeighbor[instance.nodes[nodesOrder[i]][j]].push_back(i);
		}
	}
}

std::vector<int> GrundySolver::getNodesColor() const{
	return nodesColor;
}

std::vector<int> GrundySolver::getNodesOrder() const{
	return nodesOrder;
}